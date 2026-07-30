#!/bin/bash -ue
rm -rf deeptmhmm_tmp XP_059096370.1
cp XP_059096370.1.faa input.faa

FASTA_PATH=$(readlink -f input.faa)
TMP_OUT=$(pwd)/deeptmhmm_tmp
FINAL_OUT=$(pwd)/XP_059096370.1

cd deeptmhmm
python predict.py \
    --fasta "${FASTA_PATH}" \
    --output-dir "${TMP_OUT}"

cd ..
mkdir -p "${FINAL_OUT}"
mv deeptmhmm_tmp/* "${FINAL_OUT}/"
rmdir deeptmhmm_tmp
