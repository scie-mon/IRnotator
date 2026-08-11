#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def map_auto(hardcoded_call: str) -> str:
    if hardcoded_call == "pass":
        return "accept"
    if hardcoded_call == "fail":
        return "reject"
    sys.exit(f"Unexpected hardcoded_call: {hardcoded_call!r} (expected pass|fail)")


def main() -> None:
    if len(sys.argv) != 6:
        sys.exit(
            "Usage: seed_decisions.py features.tsv hardcoded.tsv "
            "candidate_summary.tsv decision.tsv manual_review"
        )

    features_path, hardcoded_path, summary_path, decision_path = map(
        Path, sys.argv[1:5]
    )
    manual_review = sys.argv[5].lower() == "true"

    features = read_tsv(features_path)
    hardcoded = {row["protein_id"]: row for row in read_tsv(hardcoded_path)}

    required_hardcoded = {
        "protein_id", "hardcoded_call", "topology_class", "hardcoded_reason"
    }
    if hardcoded:
        missing = required_hardcoded - set(next(iter(hardcoded.values())))
        if missing:
            sys.exit(
                "hardcoded_scores.tsv missing required column(s): "
                + ", ".join(sorted(missing))
            )

    with summary_path.open("w", newline="") as sum_fh, decision_path.open(
        "w", newline=""
    ) as dec_fh:
        sum_writer = csv.writer(sum_fh, delimiter="\t")
        dec_writer = csv.writer(dec_fh, delimiter="\t")

        sum_writer.writerow([
            "protein_id",
            "feature_order",
            "n_tmhelix",
            "has_signal",
            "has_unknown_features",
            "hardcoded_call",
            "topology_class",
            "hardcoded_reason",
        ])
        dec_writer.writerow([
            "protein_id",
            "auto_decision",
            "final_decision",
            "topology_class",
            "decision_notes",
        ])

        for row in features:
            pid = row["protein_id"]
            hc = hardcoded.get(pid)
            if hc is None:
                sys.exit(f"Missing hardcoded row for protein_id={pid}")

            hardcoded_call = hc["hardcoded_call"]
            topology_class = hc["topology_class"]
            hardcoded_reason = hc.get("hardcoded_reason", "")

            sum_writer.writerow([
                pid,
                row["feature_order"],
                row["n_tmhelix"],
                row["has_signal"],
                row["has_unknown_features"],
                hardcoded_call,
                topology_class,
                hardcoded_reason,
            ])

            auto_decision = map_auto(hardcoded_call)
            final_decision = "review" if manual_review else auto_decision
            dec_writer.writerow([
                pid,
                auto_decision,
                final_decision,
                topology_class,
                "",
            ])


if __name__ == "__main__":
    main()
