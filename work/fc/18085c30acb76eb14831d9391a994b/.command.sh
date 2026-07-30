#!/bin/bash -ue
python /home/primeline/Data/IRnotator/bin/apply_manual_decisions.py \
    decision.tsv \
    review.tsv \
    decision.tsv
