process PREPARE_MANUAL_REVIEW {
    tag "prepare_manual_review"

    publishDir "${params.outdir}y", mode: 'copy'

    input:
    path deeptmhmm_dirs
    path scores_tsv
    path reviewer_html

    output:
    path "manual_review", emit: review_dir
    path "manual_review/review-manifest.json", emit: manifest
    path "manual_review/topology-reviewer.html", emit: reviewer

    script:
    """
    set -euo pipefail

    mkdir -p manual_review/plots
    cp "${reviewer_html}" manual_review/topology-reviewer.html
    cp "${scores_tsv}" manual_review/hardcoded_scores.tsv

    python3 - <<'PY'
from pathlib import Path
import csv
import json
import shutil

root = Path(".")
out_plots = Path("manual_review/plots")
out_plots.mkdir(parents=True, exist_ok=True)

seed = {}
with open("manual_review/hardcoded_scores.tsv", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    fields = list(reader.fieldnames or [])
    id_key = next((k for k in fields if k.lower() in ("seq_id", "id", "sequence_id", "protein_id")), fields[0])
    hard_key = next((k for k in fields if k.lower() in ("hardcoded_call", "hardcoded")), None)
    if hard_key is None:
        raise SystemExit(f"No 'hardcoded' column in scores TSV. Columns: {fields}")
    for row in reader:
        sid = (row.get(id_key) or "").strip()
        hard = (row.get(hard_key) or "").strip().lower()
        if hard == "pass":
            dec = "accept"
        elif hard == "fail":
            dec = "reject"
        else:
            dec = "undecided"
        if sid:
            seed[sid] = dec

items = []
for d in sorted(root.iterdir()):
    if not d.is_dir() or d.name == "manual_review":
        continue
    plot = d / "plot.png"
    if not plot.exists():
        cands = list(d.rglob("plot.png"))
        plot = cands[0] if cands else None
    if plot is None or not Path(plot).exists():
        continue
    seq_id = d.name
    shutil.copy2(plot, out_plots / f"{seq_id}.png")
    items.append({
        "id": seq_id,
        "image": f"plots/{seq_id}.png",
        "seed": seed.get(seq_id, "undecided"),
    })

if not items:
    raise SystemExit("PREPARE_MANUAL_REVIEW: no plot.png found")

Path("manual_review/review-manifest.json").write_text(
    json.dumps({"items": items}, indent=2) + "\\n"
)
n = sum(1 for i in items if i["seed"] != "undecided")
print(f"Prepared {len(items)} plots; non-undecided seeds: {n}")
PY

    cat > manual_review/README.txt <<'EOF'
Manual review package
=====================

1. From this directory, serve locally:
   python3 -m http.server 8000

2. Open:
   http://localhost:8000/topology-reviewer.html

3. Review every plot. Export review.tsv with key e.

4. Resume later with:
   --review_tsv /path/to/review.tsv -resume
EOF
    """
}
