#!/usr/bin/env python3
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Iterator

ALLOWED = {"signal", "outside", "inside", "tmhelix"}

def blocks(lines):
    buf = []
    for ln in lines:
        if ln.strip() == "//":
            if buf:
                yield buf
            buf = []
        else:
            buf.append(ln)
    if buf:
        yield buf

def normalise(ftype: str) -> str:
    f = ftype.lower()
    if f.startswith("signal"):
        return "signal"
    if f in ALLOWED:
        return f
    return f

def protein_id(block) -> str:
    for ln in block:
        if ln.startswith("#"):
            stripped = ln[1:].strip()
            if stripped:
                return stripped.split()[0]
        else:
            cols = ln.rstrip("\n").split("\t")
            if cols:
                return cols[0]
    return "unknown"

def feature_types(block):
    types = []
    for ln in block:
        if ln.startswith("#") or not ln.strip():
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 2:
            continue
        types.append(normalise(cols[1]))
    return types

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: parse_deeptmhmm_features.py TMRs.gff3 deeptmhmm_features.tsv")

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    lines = in_path.read_text().splitlines(keepends=True)
    if lines and lines[0].lstrip().startswith("##gff-version"):
        lines = lines[1:]

    with out_path.open("w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow([
            "protein_id",
            "feature_order",
            "n_tmhelix",
            "has_signal",
            "has_unknown_features"
        ])

        for block in blocks(lines):
            pid = protein_id(block)
            types = feature_types(block)
            n_tmhelix = sum(1 for t in types if t == "tmhelix")
            has_signal = "signal" in types
            has_unknown = any(t not in ALLOWED for t in types)
            writer.writerow([
                pid,
                ",".join(types),
                n_tmhelix,
                str(has_signal).lower(),
                str(has_unknown).lower()
            ])

if __name__ == "__main__":
    main()
