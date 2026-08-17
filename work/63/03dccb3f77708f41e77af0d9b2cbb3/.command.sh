#!/bin/bash -ue
set -euo pipefail

python3 merge_ir_isoforms.py \
    --gff Dmagna_LRV_IRs_braker_fixed.gff Dmagna_LRV_rnd1.gff GCA_030254905.1_UOB_LRV0_1_genomic.gff Dmagna_LRV_snap-aug.gff \
    --registry sequence_registry.tsv \
    --decisions decision.tsv \
    --overlap-fraction 0.60 \
    --mode unfiltered \
    --out-gff merged.unfiltered.gff3 \
    --out-audit merged.unfiltered.audit.tsv
