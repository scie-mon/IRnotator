#!/bin/bash -ue
set -euo pipefail

mkdir -p manual_review/plots

cp "topology-reviewer.html" manual_review/topology-reviewer.html
cp "hardcoded_scores.tsv" manual_review/hardcoded_scores.tsv
cp "sequence_registry.tsv" manual_review/sequence_registry.tsv

python3 - <<'PY'
from pathlib import Path
import csv
import json
import shutil


ROOT = Path(".")
REVIEW_DIR = Path("manual_review")
PLOTS_DIR = REVIEW_DIR / "plots"


def fail(message):
    raise SystemExit(f"PREPARE_MANUAL_REVIEW: {message}")


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            fail(f"empty TSV: {path}")
        return list(reader), reader.fieldnames


registry_rows, registry_fields = read_tsv(
    REVIEW_DIR / "sequence_registry.tsv"
)

required_registry_fields = {
    "internal_id",
    "source_id",
    "source_header",
}
missing_registry_fields = required_registry_fields - set(registry_fields)
if missing_registry_fields:
    fail(
        "sequence registry missing required column(s): "
        + ", ".join(sorted(missing_registry_fields))
    )

registry_by_internal_id = {}
for row in registry_rows:
    internal_id = (row["internal_id"] or "").strip()
    if not internal_id:
        fail("registry contains an empty internal_id")

    if internal_id in registry_by_internal_id:
        fail(f"duplicate internal_id in registry: {internal_id!r}")

    registry_by_internal_id[internal_id] = row


score_rows, score_fields = read_tsv(
    REVIEW_DIR / "hardcoded_scores.tsv"
)

required_score_fields = {
    "protein_id",
    "hardcoded_call",
    "hardcoded_reason",
}
missing_score_fields = required_score_fields - set(score_fields)
if missing_score_fields:
    fail(
        "hardcoded scores missing required column(s): "
        + ", ".join(sorted(missing_score_fields))
    )

seed_by_internal_id = {}
reason_by_internal_id = {}

for row in score_rows:
    internal_id = (row["protein_id"] or "").strip()
    hardcoded_call = (row["hardcoded_call"] or "").strip().lower()
    hardcoded_reason = (row["hardcoded_reason"] or "").strip()

    if not internal_id:
        fail("hardcoded scores contains an empty protein_id")
    if internal_id in seed_by_internal_id:
        fail(f"duplicate protein_id in hardcoded scores: {internal_id!r}")
    if internal_id not in registry_by_internal_id:
        fail(
            f"hardcoded score ID is absent from sequence registry: "
            f"{internal_id!r}"
        )

    if hardcoded_call == "pass":
        seed = "accept"
    elif hardcoded_call == "fail":
        seed = "reject"
    else:
        fail(
            f"invalid hardcoded_call for {internal_id!r}: "
            f"{hardcoded_call!r} (expected pass or fail)"
        )

    seed_by_internal_id[internal_id] = seed
    reason_by_internal_id[internal_id] = hardcoded_reason


plot_by_internal_id = {}

for directory in sorted(ROOT.iterdir()):
    if not directory.is_dir() or directory.name == REVIEW_DIR.name:
        continue

    internal_id = directory.name

    if internal_id not in registry_by_internal_id:
        fail(
            f"DeepTMHMM plot directory does not match an internal_id "
            f"in the registry: {internal_id!r}"
        )

    candidates = list(directory.rglob("plot.png"))
    if len(candidates) != 1:
        fail(
            f"expected exactly one plot.png for {internal_id!r}; "
            f"found {len(candidates)}"
        )

    if internal_id in plot_by_internal_id:
        fail(f"duplicate plot directory for internal_id: {internal_id!r}")

    plot_by_internal_id[internal_id] = candidates[0]


if not plot_by_internal_id:
    fail("no DeepTMHMM plot.png files found")

plot_ids = set(plot_by_internal_id)
score_ids = set(seed_by_internal_id)

missing_plots = score_ids - plot_ids
extra_plots = plot_ids - score_ids

if missing_plots:
    fail(
        "hardcoded-score ID(s) have no topology plot, e.g. "
        + repr(sorted(missing_plots)[0])
    )

if extra_plots:
    fail(
        "topology plot ID(s) have no hardcoded score, e.g. "
        + repr(sorted(extra_plots)[0])
    )


items = []

for internal_id in sorted(plot_by_internal_id):
    registry_row = registry_by_internal_id[internal_id]
    output_image = PLOTS_DIR / f"{internal_id}.png"

    shutil.copy2(plot_by_internal_id[internal_id], output_image)

    items.append(
        {
            "seq_id": internal_id,
            "source_id": registry_row["source_id"],
            "source_header": registry_row["source_header"],
            "image": f"plots/{internal_id}.png",
            "seed": seed_by_internal_id[internal_id],
            "reason": reason_by_internal_id[internal_id],
        }
    )

(REVIEW_DIR / "review-manifest.json").write_text(
    json.dumps({"items": items}, indent=2) + "\n"
)

accepted = sum(item["seed"] == "accept" for item in items)
rejected = sum(item["seed"] == "reject" for item in items)

print(
    f"Prepared {len(items)} plots; "
    f"accept seeds: {accepted}; reject seeds: {rejected}"
)
PY

cat > manual_review/README.txt <<'EOF'
Manual review package
=====================

Identity rules
--------------
- seq_id is the IRnotator internal ID and is the only decision key.
- source_id/source_header are displayed only for human interpretation.
- image filenames use internal IDs.
- review.tsv must retain the exported seq_id values exactly.

Instructions
------------
1. From this directory, serve locally:
   python3 -m http.server 8000

2. Open:
   http://localhost:8000/topology-reviewer.html

3. Review every plot and export review.tsv with e.

4. Resume the pipeline after review.tsv has been written here.
EOF
