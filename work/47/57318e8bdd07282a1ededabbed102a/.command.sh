#!/bin/bash -ue
python - <<'PY'
from pathlib import Path

raw = ['home', 'primeline', 'Data', 'IRnotator', 'work', '3b', '99ab72c812de68c04248cbe3b8eb88', 'deeptmhmm_summary']
paths = [Path(raw[i+1]) for i in range(0, len(raw), 2)]
out_gff = Path("merged_TMRs.gff3")

gff_files = []
for d in paths:
    gff = d / "TMRs.gff3"
    if gff.exists():
        gff_files.append(gff)

if not gff_files:
    raise SystemExit("No TMRs.gff3 files found in collected DeepTMHMM result directories")

with out_gff.open("w") as out:
    wrote_header = False
    for gff in sorted(gff_files):
        with gff.open() as fh:
            for line in fh:
                if line.startswith("#"):
                    if not wrote_header:
                        out.write(line)
                    continue
                out.write(line)
        wrote_header = True
PY

python /home/primeline/Data/IRnotator/bin/parse_deeptmhmm_features.py \
    merged_TMRs.gff3 \
    deeptmhmm_features.tsv
