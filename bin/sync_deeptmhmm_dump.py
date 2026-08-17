#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import os
import shutil
import tempfile
from pathlib import Path


def normalise_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper().rstrip("*")


def result_sequence(directory: Path) -> str | None:
    gff = [path for path in directory.glob("TMRs.gff3") if path.is_file()]
    topology = [path for path in directory.glob("predicted_topologies.3line") if path.is_file()]
    png = [path for path in directory.glob("*.png") if path.is_file()]
    if len(gff) != 1 or len(topology) != 1 or len(png) != 1:
        return None

    lines = topology[0].read_text().splitlines()
    if len(lines) < 2 or not lines[0].startswith(">"):
        return None
    return normalise_sequence(lines[1]) or None


def read_registry(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {row["sequence_sha256"]: row for row in csv.DictReader(handle, delimiter="\t")}


def write_registry(path: Path, rows: dict[str, dict[str, str]]) -> None:
    fieldnames = ["sequence_sha256", "aa_sequence", "source_sequence_id", "result_dir"]
    with tempfile.NamedTemporaryFile("w", dir=path.parent, delete=False, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows[key] for key in sorted(rows))
        temporary_path = Path(handle.name)
    os.replace(temporary_path, path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add valid DeepTMHMM result directories to a sequence-keyed local dump."
    )
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("dump_dir", type=Path)
    args = parser.parse_args()

    if not args.results_dir.is_dir():
        raise NotADirectoryError(f"DeepTMHMM results directory not found: {args.results_dir}")

    args.dump_dir.mkdir(parents=True, exist_ok=True)
    registry_path = args.dump_dir / "registry.tsv"
    lock_path = args.dump_dir / ".registry.lock"

    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        registry = read_registry(registry_path)

        for source_dir in sorted((path for path in args.results_dir.iterdir() if path.is_dir()), key=lambda path: path.name):
            sequence = result_sequence(source_dir)
            if sequence is None:
                continue

            digest = hashlib.sha256(sequence.encode()).hexdigest()
            destination = args.dump_dir / digest
            existing = registry.get(digest)

            if existing is not None:
                if existing["aa_sequence"] != sequence:
                    raise RuntimeError(f"SHA-256 collision for DeepTMHMM dump entry: {digest}")
                continue

            if destination.exists():
                existing_sequence = result_sequence(destination)
                if existing_sequence != sequence:
                    raise RuntimeError(f"Existing dump directory does not match {digest}: {destination}")
            else:
                shutil.copytree(source_dir, destination)

            registry[digest] = {
                "sequence_sha256": digest,
                "aa_sequence": sequence,
                "source_sequence_id": source_dir.name,
                "result_dir": destination.name,
            }

        write_registry(registry_path, registry)
        fcntl.flock(lock_handle, fcntl.LOCK_UN)


if __name__ == "__main__":
    main()
