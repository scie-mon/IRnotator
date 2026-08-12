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
