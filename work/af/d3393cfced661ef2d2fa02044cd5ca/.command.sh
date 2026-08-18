#!/bin/bash -ue
python /home/primeline/Data/IRnotator/bin/extract_passed_fasta.py \
    decision.tsv \
    HMM_hits.aa \
    passed.faa
