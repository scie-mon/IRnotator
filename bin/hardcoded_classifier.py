#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

RELAXED_CORE = [
    "outside", "tmhelix", "inside", "tmhelix", "outside", "tmhelix"
]
STRICT_PATTERN = [
    "signal", "outside", "tmhelix", "inside",
    "tmhelix", "outside", "tmhelix", "inside",
]


def classify(feature_order: str, has_unknown: str) -> tuple[str, str, str]:
    """Return hardcoded_call, topology_class, reason."""
    types = [x for x in feature_order.split(",") if x]

    if has_unknown == "true":
        return "fail", "fail", "unknown_feature_type"
    if not types:
        return "fail", "fail", "empty_topology"

    if types == STRICT_PATTERN:
        return "pass", "strict", "matches_strict_topology"

    core = list(types)
    if core and core[0] == "signal":
        core.pop(0)
    if core and core[-1] == "inside":
        core.pop()

    if core == RELAXED_CORE:
        return "pass", "relaxed", "matches_relaxed_topology"

    return "fail", "fail", "wrong_feature_order"


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: hardcoded_classifier.py deeptmhmm_features.tsv hardcoded_scores.tsv"
        )

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    with in_path.open() as in_fh, out_path.open("w", newline="") as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow([
            "protein_id",
            "hardcoded_call",
            "topology_class",
            "hardcoded_reason",
        ])

        for row in reader:
            call, topology_class, reason = classify(
                row["feature_order"], row["has_unknown_features"]
            )
            writer.writerow([row["protein_id"], call, topology_class, reason])


if __name__ == "__main__":
    main()
