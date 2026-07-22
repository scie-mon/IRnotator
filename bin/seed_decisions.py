#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: seed_decisions.py features.tsv hardcoded.tsv candidate_summary.tsv decision.tsv")

    features_path, hardcoded_path, summary_path, decision_path = map(Path, sys.argv[1:5])

    features = read_tsv(features_path)
    hardcoded = {row["protein_id"]: row for row in read_tsv(hardcoded_path)}

    with summary_path.open("w", newline="") as sum_fh, decision_path.open("w", newline="") as dec_fh:
        sum_writer = csv.writer(sum_fh, delimiter="\t")
        dec_writer = csv.writer(dec_fh, delimiter="\t")

        sum_writer.writerow([
            "protein_id",
            "feature_order",
            "n_tmhelix",
            "has_signal",
            "has_unknown_features",
            "hardcoded_call",
            "hardcoded_reason"
        ])

        dec_writer.writerow([
            "protein_id",
            "auto_decision",
            "final_decision",
            "decision_notes"
        ])

        for row in features:
            pid = row["protein_id"]
            hc = hardcoded.get(pid, {})
            hardcoded_call = hc.get("hardcoded_call", "na")
            hardcoded_reason = hc.get("hardcoded_reason", "na")

            sum_writer.writerow([
                pid,
                row["feature_order"],
                row["n_tmhelix"],
                row["has_signal"],
                row["has_unknown_features"],
                hardcoded_call,
                hardcoded_reason
            ])

            auto_decision = "review"
            dec_writer.writerow([pid, auto_decision, auto_decision, ""])

if __name__ == "__main__":
    main()
