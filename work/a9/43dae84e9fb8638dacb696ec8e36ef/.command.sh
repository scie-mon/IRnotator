#!/bin/bash -ue
cat > generated_results.txt <<'EOF'

EOF

    python /home/primeline/Data/IRnotator/bin/finalize_deeptmhmm_results.py \
        deeptmhmm_results.initial.tsv \
        generated_results.txt \
        deeptmhmm_results.tsv \
        deeptmhmm_results \
        unresolved_deeptmhmm.faa \
        deeptmhmm_reuse_summary.tsv

    python /home/primeline/Data/IRnotator/bin/sync_deeptmhmm_dump.py \
        deeptmhmm_results \
        /home/primeline/Data/IRnotator/deeptmhmm_dump

    if grep -q '^>' unresolved_deeptmhmm.faa; then
        echo "WARNING: DeepTMHMM results are incomplete." >&2
        echo "Unresolved sequences were written to:" >&2
        echo "  wrapper_test2/deeptmhmm/unresolved_deeptmhmm.faa" >&2
        echo "Add completed results under --deeptmhmm_salvage_paths and rerun with -resume." >&2
    else
        mv deeptmhmm_results.tsv deeptmhmm_results.complete.tsv
        mv deeptmhmm_results deeptmhmm_results.complete
    fi
