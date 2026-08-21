#!/usr/bin/env python3
"""Synchronize and sanitize a sequence-keyed DeepTMHMM result cache.

Usage:
    sync_deeptmhmm_dump.py RESULTS_DIR DUMP_DIR

RESULTS_DIR contains newly finalized DeepTMHMM result directories.  DUMP_DIR is
an on-disk cache whose directory names are SHA-256 hashes of normalized amino
acid sequences, accompanied by registry.tsv.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Iterable

REGISTRY_NAME = "registry.tsv"
LOCK_NAME = ".sync.lock"
REGISTRY_FIELDS = [
    "sequence_sha256",
    "aa_sequence",
    "source_sequence_id",
    "result_dir",
]
REQUIRED_RESULT_FILES = ("TMRs.gff3", "predicted_topologies.3line")


def normalise_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper().rstrip("*")


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def read_topology_sequence(path: Path) -> str | None:
    """Read the sole sequence from a DeepTMHMM .3line topology file."""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) != 3 or not lines[0].startswith(">"):
        return None

    sequence = normalise_sequence(lines[1])
    topology = "".join(lines[2].split())
    if not sequence or not topology or len(sequence) != len(topology):
        return None
    return sequence


def result_sequence(directory: Path) -> str | None:
    """Return the input protein for a complete single-result DeepTMHMM directory."""
    if not directory.is_dir() or directory.is_symlink():
        return None

    if any(not (directory / filename).is_file() for filename in REQUIRED_RESULT_FILES):
        return None

    pngs = [path for path in directory.glob("*.png") if path.is_file()]
    if len(pngs) != 1:
        return None

    return read_topology_sequence(directory / "predicted_topologies.3line")


def read_registry(path: Path) -> dict[str, dict[str, str]]:
    """Read valid registry rows, ignoring malformed or internally inconsistent rows."""
    if not path.exists():
        return {}
    if not path.is_file():
        raise RuntimeError(f"Registry path is not a regular file: {path}")

    registry: dict[str, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REGISTRY_FIELDS:
            raise RuntimeError(
                f"Unexpected registry header in {path}; expected: {', '.join(REGISTRY_FIELDS)}"
            )

        for line_number, row in enumerate(reader, start=2):
            if row is None or any(row.get(field) is None for field in REGISTRY_FIELDS):
                print(f"WARNING: ignoring malformed registry row {line_number}", file=sys.stderr)
                continue

            sequence = normalise_sequence(row["aa_sequence"])
            digest = sequence_sha256(sequence) if sequence else ""
            result_dir = row["result_dir"].strip()
            if (
                not sequence
                or row["sequence_sha256"] != digest
                or not result_dir
                or Path(result_dir).name != result_dir
            ):
                print(f"WARNING: ignoring invalid registry row {line_number}", file=sys.stderr)
                continue

            candidate = {
                "sequence_sha256": digest,
                "aa_sequence": sequence,
                "source_sequence_id": row["source_sequence_id"].strip(),
                "result_dir": result_dir,
            }
            existing = registry.get(digest)
            if existing is None or candidate["result_dir"] < existing["result_dir"]:
                registry[digest] = candidate
    return registry


def write_registry(path: Path, rows: dict[str, dict[str, str]]) -> None:
    """Atomically replace the registry in deterministic hash order."""
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
    )
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=REGISTRY_FIELDS, delimiter="\t")
            writer.writeheader()
            for digest in sorted(rows):
                writer.writerow(rows[digest])
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def remove_entry(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def copy_result(source: Path, destination: Path) -> None:
    """Copy a validated result without following a pre-existing destination."""
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to overwrite cache directory: {destination}")

    temporary = destination.parent / f".{destination.name}.copying"
    remove_entry(temporary)
    try:
        shutil.copytree(source, temporary, symlinks=True)
        os.replace(temporary, destination)
    except BaseException:
        remove_entry(temporary)
        raise


def valid_dump_entries(dump_dir: Path) -> Iterable[tuple[Path, str]]:
    for entry in sorted(dump_dir.iterdir(), key=lambda item: item.name):
        if entry.name in {REGISTRY_NAME, LOCK_NAME} or entry.name.startswith(f".{REGISTRY_NAME}."):
            continue
        sequence = result_sequence(entry)
        if sequence is not None:
            yield entry, sequence


def choose_canonical(
    digest: str,
    sequence: str,
    registry_row: dict[str, str] | None,
    dump_candidates: list[Path],
    source_candidates: list[Path],
) -> tuple[dict[str, str], Path | None]:
    """Select the stable canonical cache result and an optional source to import."""
    dump_by_name = {path.name: path for path in dump_candidates}
    source_by_name = {path.name: path for path in source_candidates}

    if registry_row is not None and registry_row["result_dir"] in dump_by_name:
        return registry_row, None

    if dump_candidates:
        chosen = min(dump_candidates, key=lambda path: path.name)
        return {
            "sequence_sha256": digest,
            "aa_sequence": sequence,
            "source_sequence_id": chosen.name,
            "result_dir": chosen.name,
        }, None

    destination = dump_by_name.get(digest)
    if destination is not None:
        raise RuntimeError(f"Unexpected pre-existing cache destination: {destination}")

    if not source_candidates:
        raise RuntimeError("No candidate result available for canonicalization")

    chosen = min(source_candidates, key=lambda path: path.name)
    return {
        "sequence_sha256": digest,
        "aa_sequence": sequence,
        "source_sequence_id": chosen.name,
        "result_dir": digest,
    }, chosen


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Synchronize and sanitize a sequence-keyed DeepTMHMM local dump."
    )
    parser.add_argument("results_dir", type=Path, help="Finalized DeepTMHMM result directory")
    parser.add_argument("dump_dir", type=Path, help="Persistent sequence-keyed DeepTMHMM cache")
    args = parser.parse_args()

    if not args.results_dir.is_dir():
        raise NotADirectoryError(f"DeepTMHMM results directory not found: {args.results_dir}")

    args.dump_dir.mkdir(parents=True, exist_ok=True)
    if not args.dump_dir.is_dir():
        raise NotADirectoryError(f"DeepTMHMM dump path is not a directory: {args.dump_dir}")

    registry_path = args.dump_dir / REGISTRY_NAME
    lock_path = args.dump_dir / LOCK_NAME

    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        try:
            registry = read_registry(registry_path)

            dump_by_hash: dict[str, list[Path]] = {}
            sequence_by_hash: dict[str, str] = {}
            for directory, sequence in valid_dump_entries(args.dump_dir):
                digest = sequence_sha256(sequence)
                dump_by_hash.setdefault(digest, []).append(directory)
                sequence_by_hash[digest] = sequence

            source_by_hash: dict[str, list[Path]] = {}
            for directory in sorted(args.results_dir.iterdir(), key=lambda item: item.name):
                sequence = result_sequence(directory)
                if sequence is None:
                    continue
                digest = sequence_sha256(sequence)
                source_by_hash.setdefault(digest, []).append(directory)
                sequence_by_hash[digest] = sequence

            reconciled: dict[str, dict[str, str]] = {}
            for digest in sorted(set(registry) | set(dump_by_hash) | set(source_by_hash)):
                sequence = sequence_by_hash.get(digest)
                if sequence is None:
                    # A stale registry row without a valid cache result or source result.
                    continue

                row, import_source = choose_canonical(
                    digest,
                    sequence,
                    registry.get(digest),
                    dump_by_hash.get(digest, []),
                    source_by_hash.get(digest, []),
                )
                if import_source is not None:
                    copy_result(import_source, args.dump_dir / row["result_dir"])
                reconciled[digest] = row

            retained = {row["result_dir"] for row in reconciled.values()}
            for entry in list(args.dump_dir.iterdir()):
                if entry.name in {REGISTRY_NAME, LOCK_NAME}:
                    continue
                if entry.name.startswith(f".{REGISTRY_NAME}."):
                    remove_entry(entry)
                    continue
                if entry.name not in retained:
                    remove_entry(entry)

            # Verify post-pruning invariants before publishing the new registry.
            for digest, row in reconciled.items():
                directory = args.dump_dir / row["result_dir"]
                sequence = result_sequence(directory)
                if sequence is None or sequence_sha256(sequence) != digest:
                    raise RuntimeError(f"Canonical dump result failed validation: {directory}")

            write_registry(registry_path, reconciled)
        finally:
            fcntl.flock(lock_handle, fcntl.LOCK_UN)


if __name__ == "__main__":
    main()
