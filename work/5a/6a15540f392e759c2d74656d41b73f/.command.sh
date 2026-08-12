#!/bin/bash -ue
cat > salvage_roots.txt <<'EOF'
/home/primeline/Data/maker/metagenome/12_Dmagna_LRV
EOF

    python /home/primeline/Data/IRnotator/bin/initialize_deeptmhmm_results.py \
        HMM_hits.aa \
        salvage_roots.txt \
        deeptmhmm_results.initial.tsv \
        unmatched \
        deeptmhmm_initial_summary.tsv
