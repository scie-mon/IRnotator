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
