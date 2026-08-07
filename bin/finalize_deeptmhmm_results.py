#!/usr/bin/env python3
from __future__ import annotations
import csv
import shutil
import sys
from pathlib import Path

def normalise_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper().rstrip("*")

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
    if len(sys.argv) != 7:
        sys.exit("Usage: finalize_deeptmhmm_results.py INITIAL.tsv GENERATED.txt FINAL.tsv RESULTS_DIR UNRESOLVED.faa SUMMARY.tsv")
    initial_path, generated_list, final_path, results_dir, unresolved_path, summary_path = map(Path, sys.argv[1:])
    generated_roots = [Path(p) for p in generated_list.read_text().splitlines() if p]
    with initial_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    missing_sequences = {row["aa_sequence"] for row in rows if not row["result_dir"]}
    generated = {}
    for directory in generated_roots:
        sequence = candidate_sequence(directory)
        if sequence in missing_sequences and sequence not in generated:
            generated[sequence] = directory.resolve()
    for row in rows:
        if not row["result_dir"] and row["aa_sequence"] in generated:
            result = generated[row["aa_sequence"]]
            row.update(result_dir=row["sequence_id"], source="generated", candidate_path=str(result))
        elif row["result_dir"]:
            row["result_dir"] = row["sequence_id"]
    results_dir.mkdir(parents=True, exist_ok=True)
    with unresolved_path.open("w") as unresolved:
        for row in rows:
            if not row["result_dir"]:
                unresolved.write(f">{row['sequence_id']}\n{row['aa_sequence']}\n")
                continue
            destination = results_dir / row["sequence_id"]
            if destination.exists():
                raise RuntimeError(f"Duplicate output result directory: {destination}")
            shutil.copytree(Path(row["candidate_path"]), destination)
    with final_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts = {
        "hmm_candidates": len(rows),
        "salvaged": sum(row["source"] == "salvaged" for row in rows),
        "generated": sum(row["source"] == "generated" for row in rows),
        "unresolved": sum(not row["result_dir"] for row in rows),
    }
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "count"])
        writer.writerows(counts.items())

if __name__ == "__main__":
    main()
