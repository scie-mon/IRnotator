#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path

ALLOWED_CORE = ["outside", "tmhelix", "inside", "tmhelix", "outside", "tmhelix"]

def classify(feature_order: str, has_unknown: str):
    types = [x for x in feature_order.split(",") if x]

    if has_unknown == "true":
        return "fail", "unknown_feature_type"

    if not types:
        return "fail", "empty_topology"

    idx = 1 if types[0] == "signal" else 0
    core = types[idx:]

    if len(core) not in (6, 7):
        return "fail", "wrong_feature_count"

    if core[:6] != ALLOWED_CORE:
        return "fail", "wrong_feature_order"

    if len(core) == 7 and core[6] != "inside":
        return "fail", "invalid_terminal_feature"

    return "pass", "matches_strict_topology"

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: hardcoded_classifier.py deeptmhmm_features.tsv hardcoded_scores.tsv")

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    with in_path.open() as in_fh, out_path.open("w", newline="") as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["protein_id", "hardcoded_call", "hardcoded_reason"])

        for row in reader:
            call, reason = classify(row["feature_order"], row["has_unknown_features"])
            writer.writerow([row["protein_id"], call, reason])

if __name__ == "__main__":
    main()
