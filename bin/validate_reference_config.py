#!/usr/bin/env python3
"""Emit non-fatal IRnotator reference-configuration diagnostics.

The script writes WARNING:/NOTE: records to stdout and always exits successfully
for configuration findings. It is suitable for direct controller invocation or
for execution from a lightweight Nextflow audit process.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

CANONICAL_HMMS = ("PF00060.hmm", "PF10613.hmm")
MANIFEST_HEADER = ("filename", "sha256")


def report(level: str, message: str) -> None:
    print(f"{level}: {message}")


def resolve_path(value: str) -> Path:
    return Path(value).expanduser().resolve(strict=False)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: Path) -> dict[str, str] | None:
    if not path.is_file():
        report("WARNING", f"Canonical HMM checksum manifest is missing: {path}")
        return None

    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != MANIFEST_HEADER:
                report("WARNING", f"Canonical HMM checksum manifest has an invalid header: {path}")
                return None

            checksums: dict[str, str] = {}
            for line_number, row in enumerate(reader, start=2):
                filename = (row.get("filename") or "").strip()
                checksum = (row.get("sha256") or "").strip().lower()
                if not filename or len(checksum) != 64 or any(char not in "0123456789abcdef" for char in checksum):
                    report("WARNING", f"Ignoring malformed manifest row {line_number}: {path}")
                    continue
                checksums[filename] = checksum
    except OSError as error:
        report("WARNING", f"Could not read canonical HMM checksum manifest {path}: {error}")
        return None

    missing = [filename for filename in CANONICAL_HMMS if filename not in checksums]
    if missing:
        report("WARNING", f"Canonical HMM checksum manifest lacks: {', '.join(missing)}")
        return None

    return {filename: checksums[filename] for filename in CANONICAL_HMMS}


def hmm_checksums(directory: Path) -> dict[Path, str]:
    if not directory.is_dir():
        return {}

    checksums: dict[Path, str] = {}
    for path in sorted(directory.iterdir(), key=lambda entry: entry.name):
        if path.is_file() and path.suffix == ".hmm":
            try:
                checksums[path] = sha256(path)
            except OSError as error:
                report("WARNING", f"Could not checksum loaded HMM {path}: {error}")
    return checksums


def audit_loaded_hmms(hmm_dir: Path, expected: dict[str, str] | None) -> None:
    if not hmm_dir.is_dir():
        report("WARNING", f"Configured HMM directory does not exist: {hmm_dir}")
        return

    loaded = hmm_checksums(hmm_dir)
    if expected is None:
        return

    loaded_hashes = set(loaded.values())
    expected_hashes = set(expected.values())
    for filename, checksum in expected.items():
        if checksum not in loaded_hashes:
            report("WARNING", f"Canonical {filename} is absent from loaded HMMs: {hmm_dir}")

    noncanonical_hashes = loaded_hashes - expected_hashes
    if noncanonical_hashes:
        report("NOTE", f"Non-canonical HMM content is in use from: {hmm_dir}")


def audit_canonical_directory(canonical_dir: Path, expected: dict[str, str] | None) -> None:
    expected_names = set(CANONICAL_HMMS)
    entries = set()
    if canonical_dir.is_dir():
        try:
            entries = {path.name for path in canonical_dir.iterdir()}
        except OSError as error:
            report("WARNING", f"Could not inspect canonical HMM directory {canonical_dir}: {error}")
            return

    layout_valid = entries == expected_names
    checksums_valid = expected is not None and layout_valid
    if checksums_valid:
        for filename, checksum in expected.items():
            path = canonical_dir / filename
            try:
                if not path.is_file() or sha256(path) != checksum:
                    checksums_valid = False
                    break
            except OSError:
                checksums_valid = False
                break

    if not layout_valid or not checksums_valid:
        report(
            "WARNING",
            f"Canonical HMM directory is incomplete, modified, or contains unexpected entries: {canonical_dir}",
        )


def audit_salvage_paths(canonical_dump: Path, salvage_paths: list[Path]) -> None:
    if canonical_dump not in salvage_paths:
        report(
            "NOTE",
            f"Default DeepTMHMM dump is not configured as a salvage path: {canonical_dump}",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Emit non-fatal diagnostics for IRnotator HMM and DeepTMHMM cache configuration."
    )
    parser.add_argument("--project-dir", required=True, type=resolve_path)
    parser.add_argument("--hmm-dir", required=True, type=resolve_path)
    parser.add_argument("--salvage-path", action="append", default=[], type=resolve_path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_dir: Path = args.project_dir
    canonical_hmm_dir = (project_dir / "hmms").resolve(strict=False)
    manifest = project_dir / "assets" / "canonical_hmms.sha256.tsv"
    canonical_dump = (project_dir / "deeptmhmm_dump").resolve(strict=False)

    expected = read_manifest(manifest)
    audit_loaded_hmms(args.hmm_dir, expected)
    if args.hmm_dir != canonical_hmm_dir:
        report("NOTE", f"Non-default HMM directory configured: {args.hmm_dir}")
    audit_canonical_directory(canonical_hmm_dir, expected)
    audit_salvage_paths(canonical_dump, args.salvage_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        report("WARNING", f"Reference-configuration audit could not complete: {error}")
        sys.exit(0)
