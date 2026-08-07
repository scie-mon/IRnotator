#!/usr/bin/env python3
from __future__ import annotations
import csv
import sys
from pathlib import Path

def normalise_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper().rstrip("*")

def read_fasta(path: Path):
    record_id, sequence = None, []
    for raw in path.open():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if record_id is not None:
                yield record_id, "".join(sequence)
            record_id, sequence = line[1:].split()[0], []
        else:
            if record_id is None:
                raise ValueError(f"Sequence data precedes a FASTA header in {path}")
            sequence.append(line)
    if record_id is not None:
        yield record_id, "".join(sequence)

def candidate_sequence(directory: Path) -> str | None:
    gff = [p for p in directory.glob("TMRs.gff3") if p.is_file()]
    topology = [p for p in directory.glob("predicted_topologies.3line") if p.is_file()]
    png = [p for p in directory.glob("*.png") if p.is_file()]
    if len(gff) != 1 or len(topology) != 1 or len(png) != 1:
        return None
    lines = topology[0].read_text().splitlines()
    if len(lines) < 2 or not lines[0].startswith(">"):
        return None
    return normalise_sequence(lines[1]) or None

def main() -> None:
    if len(sys.argv) != 6:
        sys.exit("Usage: initialize_deeptmhmm_results.py INPUT.faa ROOTS.txt MANIFEST.tsv UNMATCHED_DIR SUMMARY.tsv")
    fasta_path, roots_path, manifest_path, unmatched_dir, summary_path = map(Path, sys.argv[1:])
    records = [(seq_id, normalise_sequence(seq)) for seq_id, seq in read_fasta(fasta_path)]
    if not records:
        raise RuntimeError("No HMM-passing sequences were supplied")
    wanted, matches = {seq for _, seq in records}, {}
    roots = [Path(x.strip()) for x in roots_path.read_text().splitlines() if x.strip()]
    for root in roots:
        if not root.is_dir():
            raise NotADirectoryError(f"DeepTMHMM salvage root is not a directory: {root}")
        directories = [root, *sorted((p for p in root.rglob("*") if p.is_dir()), key=str)]
        for directory in directories:
            if len(matches) == len(wanted):
                break
            sequence = candidate_sequence(directory)
            if sequence in wanted and sequence not in matches:
                matches[sequence] = directory.resolve()
        if len(matches) == len(wanted):
            break
    unmatched_dir.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence_id", "aa_sequence", "result_dir", "source", "candidate_path"], delimiter="\t")
        writer.writeheader()
        for sequence_id, sequence in records:
            candidate = matches.get(sequence)
            if candidate is None:
                (unmatched_dir / f"{sequence_id}.faa").write_text(f">{sequence_id}\n{sequence}\n")
                writer.writerow(dict(sequence_id=sequence_id, aa_sequence=sequence, result_dir="", source="", candidate_path=""))
            else:
                writer.writerow(dict(sequence_id=sequence_id, aa_sequence=sequence, result_dir=str(candidate), source="salvaged", candidate_path=str(candidate)))
    salvaged = sum(1 for _, sequence in records if sequence in matches)
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "count"])
        writer.writerows([
            ["hmm_candidates", len(records)],
            ["salvaged", salvaged],
            ["queued_for_prediction", len(records) - salvaged],
        ])

if __name__ == "__main__":
    main()
