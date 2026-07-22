#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

def read_ids(path: Path) -> set[str]:
    with path.open() as fh:
        return {line.strip() for line in fh if line.strip()}

def write_filtered(ids: set[str], fasta_in: Path, fasta_out: Path) -> None:
    with fasta_in.open() as fin, fasta_out.open("w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                keep = seq_id in ids
            if keep:
                fout.write(line)

def main(argv: list[str]) -> None:
    if len(argv) != 4:
        sys.exit("Usage: python extract_fasta_by_id.py ids.txt input.faa output.faa")

    ids_path, fasta_in_path, fasta_out_path = map(Path, argv[1:])
    ids = read_ids(ids_path)
    write_filtered(ids, fasta_in_path, fasta_out_path)

if __name__ == "__main__":
    main(sys.argv)
