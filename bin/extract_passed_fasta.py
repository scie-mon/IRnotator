#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path


def load_accept_ids(decision_path: Path) -> list[str]:
    ids: list[str] = []
    with decision_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "protein_id" not in (reader.fieldnames or []):
            sys.exit(f"decision.tsv missing protein_id column: {reader.fieldnames}")
        if "final_decision" not in (reader.fieldnames or []):
            sys.exit(f"decision.tsv missing final_decision column: {reader.fieldnames}")
        for row in reader:
            if row["final_decision"].strip() == "accept":
                pid = row["protein_id"].strip()
                if pid:
                    ids.append(pid)
    return ids


def iter_fasta(path: Path):
    header = None
    seq_chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0] if line[1:] else ""
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: extract_passed_fasta.py decision.tsv hmm_hits.faa passed.faa"
        )

    decision_path, fasta_path, out_path = map(Path, sys.argv[1:4])
    accept_ids = load_accept_ids(decision_path)
    accept_set = set(accept_ids)

    written = 0
    seen_in_fasta: set[str] = set()

    with out_path.open("w") as out_fh:
        for seq_id, seq in iter_fasta(fasta_path):
            if seq_id in accept_set:
                out_fh.write(f">{seq_id}\n")
                for i in range(0, len(seq), 60):
                    out_fh.write(seq[i : i + 60] + "\n")
                written += 1
                seen_in_fasta.add(seq_id)

    missing = accept_set - seen_in_fasta
    if missing:
        example = next(iter(missing))
        sys.exit(
            f"{len(missing)} accepted id(s) not found in FASTA, e.g. {example}"
        )

    print(f"EXTRACT_PASSED_FASTA: wrote {written} sequence(s) to {out_path.name}")


if __name__ == "__main__":
    main()
