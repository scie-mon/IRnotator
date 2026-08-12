#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit("Usage: apply_manual_decisions.py seed_decision.tsv review.tsv decision.tsv")

    seed_path, review_path, out_path = map(Path, sys.argv[1:4])
    seed_rows = read_tsv(seed_path)
    review_rows = read_tsv(review_path)

    if seed_rows:
        required = {"protein_id", "auto_decision", "topology_class"}
        missing = required - set(seed_rows[0])
        if missing:
            sys.exit("seed decision TSV missing required column(s): " + ", ".join(sorted(missing)))

    review_by_id: dict[str, str] = {}
    for row in review_rows:
        sequence_id = row.get("seq_id") or row.get("protein_id")
        if not sequence_id:
            sys.exit(f"Review row missing seq_id: {row}")
        decision = row.get("decision", "").strip()
        if decision not in ("accept", "reject"):
            sys.exit(f"Invalid decision for {sequence_id!r}: {decision!r} (expected accept|reject)")
        if sequence_id in review_by_id:
            sys.exit(f"Duplicate seq_id in review.tsv: {sequence_id}")
        review_by_id[sequence_id] = decision

    seed_ids = [row["protein_id"] for row in seed_rows]
    seed_set = set(seed_ids)
    review_set = set(review_by_id)
    extra = review_set - seed_set
    if extra:
        sys.exit(f"review.tsv has {len(extra)} unknown id(s), e.g. {next(iter(extra))}")

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "protein_id",
            "auto_decision",
            "final_decision",
            "topology_class",
            "decision_notes",
        ])
        for row in seed_rows:
            protein_id = row["protein_id"]
            writer.writerow([
                protein_id,
                row["auto_decision"],
                review_by_id.get(protein_id, "reject"),
                row["topology_class"],
                row.get("decision_notes", ""),
            ])


if __name__ == "__main__":
    main()
