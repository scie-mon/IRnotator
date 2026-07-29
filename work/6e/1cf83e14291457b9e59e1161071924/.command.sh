#!/bin/bash -ue
rm -rf deeptmhmm_tmp LOC116915405
cp LOC116915405.faa input.faa

FASTA_PATH=$(readlink -f input.faa)
TMP_OUT=$(pwd)/deeptmhmm_tmp
FINAL_OUT=$(pwd)/LOC116915405

cd deeptmhmm
python predict.py \
    --fasta "${FASTA_PATH}" \
    --output-dir "${TMP_OUT}"

cd ..
mkdir -p "${FINAL_OUT}"
mv deeptmhmm_tmp/* "${FINAL_OUT}/"
rmdir deeptmhmm_tmp
