#!/bin/bash -ue
set -euo pipefail

    mkdir -p manual_review/plots manual_review/genes
    cp topology-reviewer.html manual_review/topology-reviewer.html
    cp hardcoded_scores.tsv manual_review/hardcoded_scores.tsv
    cp sequence_registry.tsv manual_review/sequence_registry.tsv
    cp deeptmhmm_features.tsv manual_review/deeptmhmm_features.tsv
    cp decision.tsv manual_review/effective_decisions.tsv
    cp merged.rev1.gff3 manual_review/reviewer_context.gff3

    python3 - <<'PY'
from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

REVIEW_DIR = Path("manual_review")
VALID_DECISIONS = {"accept", "reject", "review"}


def fail(message):
    raise SystemExit(f"PREPARE_MANUAL_REVIEW: {message}")


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            fail(f"empty TSV: {path}")
        return list(reader)


def attrs(text):
    return {
        key: value
        for item in text.rstrip(";").split(";")
        if "=" in item
        for key, value in [item.split("=", 1)]
    }


def parse_context(path):
    items = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                fail(f"{path}:{line_number}: expected 9-column GFF3")
            contig, _, kind, start, end, _, strand, _, raw_attrs = fields
            if kind not in {"mRNA", "transcript"}:
                continue
            parsed = attrs(raw_attrs)
            internal_id = parsed.get("internal_id", "").strip()
            transcript_id = parsed.get("ID", "").strip()
            parents = [parent for parent in parsed.get("Parent", "").split(",") if parent]
            if not internal_id or not transcript_id or len(parents) != 1:
                fail(f"{path}:{line_number}: mRNA requires internal_id, ID, and one Parent")
            if internal_id in items:
                fail(f"duplicate internal_id in reviewer context GFF: {internal_id}")
            items[internal_id] = {
                "seq_id": internal_id,
                "gene_id": parents[0],
                "transcript_id": transcript_id,
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
            }
    if not items:
        fail(f"no mRNA records found in reviewer context GFF: {path}")
    return items


def normalize_decision(decision):
    return "undecided" if decision == "review" else decision


registry_by_id = {}
for row in read_tsv(REVIEW_DIR / "sequence_registry.tsv"):
    internal_id = row.get("internal_id", "").strip()
    if not internal_id or internal_id in registry_by_id:
        fail(f"invalid or duplicate internal_id in sequence registry: {internal_id!r}")
    registry_by_id[internal_id] = row

reason_by_id = {}
for row in read_tsv(REVIEW_DIR / "hardcoded_scores.tsv"):
    internal_id = row.get("protein_id", "").strip()
    if internal_id:
        reason_by_id[internal_id] = row.get("hardcoded_reason", "").strip()

decision_by_id = {}
for row in read_tsv(REVIEW_DIR / "effective_decisions.tsv"):
    internal_id = row.get("protein_id", "").strip()
    decision = row.get("final_decision", "").strip().lower()
    if not internal_id:
        fail("empty protein_id in effective decisions TSV")
    if decision not in VALID_DECISIONS:
        fail(f"invalid final_decision for {internal_id}: {decision!r}")
    if internal_id in decision_by_id:
        fail(f"duplicate protein_id in effective decisions TSV: {internal_id}")
    decision_by_id[internal_id] = decision

items = []
for internal_id, item in parse_context(REVIEW_DIR / "reviewer_context.gff3").items():
    if internal_id not in registry_by_id:
        fail(f"reviewer-GFF internal_id is absent from registry: {internal_id}")
    if internal_id not in decision_by_id:
        fail(f"reviewer-GFF internal_id is absent from effective decisions: {internal_id}")

    plots = list(Path(internal_id).rglob("*.png"))
    if len(plots) != 1:
        fail(f"expected one DeepTMHMM PNG for {internal_id}; found {len(plots)}")

    plot_out = REVIEW_DIR / "plots" / f"{internal_id}.png"
    gene_out = REVIEW_DIR / "genes" / f"{internal_id}.svg"
    shutil.copy2(plots[0], plot_out)
    subprocess.run([
        "python3", "render_gene_cds.py", "merged.rev1.gff3", str(gene_out),
        "--features", "deeptmhmm_features.tsv",
        "--decisions", "decision.tsv",
        "--focal-internal-id", internal_id,
    ], check=True)

    registry = registry_by_id[internal_id]
    item["seed"] = normalize_decision(decision_by_id[internal_id])
    item.update({
        "source_id": registry.get("source_id", ""),
        "source_header": registry.get("source_header", ""),
        "reason": reason_by_id.get(internal_id, ""),
        "topology_image": f"plots/{internal_id}.png",
        "gene_image": f"genes/{internal_id}.svg",
    })
    items.append(item)

items.sort(key=lambda item: (
    item["contig"], item["start"], item["gene_id"], item["end"], item["transcript_id"],
))
for order, item in enumerate(items, 1):
    item["review_order"] = order

(REVIEW_DIR / "review-manifest.json").write_text(json.dumps({
    "reviewer_revision": 1,
    "items": items,
}, indent=2) + "\n")
PY

    if [[ 2 -gt 0 ]]; then
        cp review.tsv manual_review/review.applied.rev2.tsv
        rm -f "test_out/manual_review/review.tsv"
    fi

    cat > manual_review/rerun_after_review.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

# Remove only cached tasks that consume review state or generate review assets.
# Usage: ./rerun_after_review.sh [WORK_DIR]
work_dir="${1:-work}"
if [[ ! -d "${work_dir}" ]]; then
    echo "Work directory not found: ${work_dir}" >&2
    exit 1
fi

find "${work_dir}" -type f  -name .command.sh -o -name .command.run  -print0 |
while IFS= read -r -d '' task_file; do
    if grep -qE 'merge_ir_isoforms.py|render_gene_cds.py|reviewer_context.gff3|review.tsv' "${task_file}"; then
        task_dir="$(dirname "${task_file}")"
        echo "Removing cached review-dependent task: ${task_dir}"
        rm -rf "${task_dir}"
    fi
done
EOF
    chmod +x manual_review/rerun_after_review.sh

    cat > manual_review/README.txt <<'EOF'
Manual review package
=====================

The manifest is ordered by merged gene model, then isoform.

1. Serve this directory: python3 -m http.server 8000
2. Open topology-reviewer.html.
3. Export review.tsv into this directory.
4. Run ./rerun_after_review.sh /path/to/IRnotator/work
5. Re-run the same Nextflow command with -resume.

seq_id is the IRnotator internal ID and is the only decision key.
EOF
