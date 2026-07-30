#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_tsv(path: Path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: apply_manual_decisions.py seed_decision.tsv review.tsv decision.tsv"
        )

    seed_path, review_path, out_path = map(Path, sys.argv[1:4])

    seed_rows = read_tsv(seed_path)
    review_rows = read_tsv(review_path)

    review_by_id = {}
    for row in review_rows:
        sid = row.get("seq_id") or row.get("protein_id")
        if not sid:
            sys.exit(f"Review row missing seq_id: {row}")
        decision = row.get("decision", "").strip()
        if decision not in ("accept", "reject"):
            sys.exit(
                f"Invalid decision for {sid!r}: {decision!r} (expected accept|reject)"
            )
        if sid in review_by_id:
            sys.exit(f"Duplicate seq_id in review.tsv: {sid}")
        review_by_id[sid] = decision

    seed_ids = [row["protein_id"] for row in seed_rows]
    seed_set = set(seed_ids)
    review_set = set(review_by_id)

    missing = seed_set - review_set
    extra = review_set - seed_set
    if missing:
        sys.exit(f"review.tsv missing {len(missing)} seed id(s), e.g. {next(iter(missing))}")
    if extra:
        sys.exit(f"review.tsv has {len(extra)} unknown id(s), e.g. {next(iter(extra))}")

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["protein_id", "auto_decision", "final_decision", "decision_notes"]
        )
        for row in seed_rows:
            pid = row["protein_id"]
            writer.writerow(
                [
                    pid,
                    row["auto_decision"],
                    review_by_id[pid],
                    "",
                ]
            )


if __name__ == "__main__":
    main()
