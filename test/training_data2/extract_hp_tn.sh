#!/usr/bin/env bash
set -euo pipefail

src_root="/home/primeline/Data/maker/metagenome"

png_target="./negative_png_target"
csv_target="./negative_csv_target"
line_target="./negative_3line_target"

mkdir -p "$png_target" "$csv_target" "$line_target"

i=1

find "$src_root" -type f -name "*.png" -path "*/HP_TN/*" -print0 |
while IFS= read -r -d '' png; do
    png_name=$(basename "$png")
    png_base="${png_name%.png}"
    png_parent=$(dirname "$png")

    results_dir="$(realpath "$png_parent/../../results")"
    probs_csv="$results_dir/${png_base}/${png_base}_probs.csv"
    topo_3line="$results_dir/${png_base}/predicted_topologies.3line"

    # all-or-none rule
    if [[ ! -f "$probs_csv" || ! -f "$topo_3line" ]]; then
        echo "Skipping incomplete set:" >&2
        echo "  PNG   : $png" >&2
        echo "  CSV   : $probs_csv" >&2
        echo "  3LINE : $topo_3line" >&2
        continue
    fi

    new_base=$(printf "%04d_%s" "$i" "$png_base")

    cp -- "$png"        "$png_target/${new_base}.png"
    cp -- "$probs_csv"  "$csv_target/${new_base}.csv"
    cp -- "$topo_3line" "$line_target/${new_base}.3line"

    ((i++))
done