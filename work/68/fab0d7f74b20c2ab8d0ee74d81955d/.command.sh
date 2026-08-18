#!/bin/bash -ue
mkdir -p deeptmhmm_summary

python - <<'PY'
import csv
from pathlib import Path

manifest = Path('deeptmhmm_results.complete.tsv')
results_root = Path('deeptmhmm_results.complete')
output = Path('deeptmhmm_summary/TMRs.gff3')

with manifest.open(newline='') as handle:
    rows = list(csv.DictReader(handle, delimiter='\t'))

if not rows:
    raise RuntimeError('No DeepTMHMM result rows to merge')

wrote_version = False
with output.open('w') as out:
    for index, row in enumerate(rows):
        sequence_id = row['sequence_id']
        result_dir = row['result_dir']
        if not result_dir:
            raise RuntimeError(f'Missing DeepTMHMM result for {sequence_id}')

        gff = results_root / result_dir / 'TMRs.gff3'
        if not gff.is_file():
            raise FileNotFoundError(f'Missing TMRs.gff3 for {sequence_id}: {gff}')

        if index:
            out.write('//\n')

        wrote_identity = False
        for line in gff.open():
            if line.startswith('##gff-version'):
                if not wrote_version:
                    out.write(line)
                    wrote_version = True
                continue
            if line.startswith('#') and not wrote_identity:
                out.write(f'# {sequence_id}\n')
                wrote_identity = True
                continue
            if not wrote_identity and line.strip():
                out.write(f'# {sequence_id}\n')
                wrote_identity = True
            out.write(line)

        if not wrote_identity:
            out.write(f'# {sequence_id}\n')
PY
