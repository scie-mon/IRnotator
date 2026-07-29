#!/bin/bash -ue
python - <<'PY'
from pathlib import Path
from Bio import SeqIO
import re

inp = Path("HMM_hits.aa")
for rec in SeqIO.parse(str(inp), "fasta"):
    seq_id = re.sub(r'[^A-Za-z0-9_.-]+', '_', rec.id)
    out = Path(f"{seq_id}.faa")
    SeqIO.write(rec, str(out), "fasta")
PY
