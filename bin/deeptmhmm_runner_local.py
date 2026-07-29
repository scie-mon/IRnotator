#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path


def iter_fasta(path: Path):
    header = None
    seq_lines = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, ''.join(seq_lines)


def safe_name(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', name)


def find_predict_script(deeptmhmm_dir: Path) -> Path:
    base = deeptmhmm_dir.resolve()
    candidates = [
        base / 'predict.py',
        base / 'deeptmhmm' / 'predict.py',
    ]
    for p in candidates:
        if p.is_file():
            return p.resolve()
    raise FileNotFoundError(f'Could not find predict.py under {base}')


def debug_package_layout(deeptmhmm_dir: Path, predict_py: Path):
    base = deeptmhmm_dir.resolve()
    print(f'[DEBUG] deeptmhmm_dir={base}', flush=True)
    print(f'[DEBUG] predict_py={predict_py}', flush=True)
    interesting = [
        base / 'esm_model_args.pt',
        base / 'deeptmhmm' / 'esm_model_args.pt',
        predict_py.parent / 'esm_model_args.pt',
        base / 'utils.py',
        base / 'deeptmhmm' / 'utils.py',
        predict_py.parent / 'utils.py',
    ]
    for p in interesting:
        print(f'[DEBUG] exists {p}: {p.exists()}', flush=True)
    try:
        top = sorted(base.iterdir())
        print('[DEBUG] top-level contents:', flush=True)
        for p in top[:50]:
            print(f'[DEBUG]   {p.name}', flush=True)
    except Exception as e:
        print(f'[DEBUG] failed listing {base}: {e}', flush=True)


def run_deeptmhmm_local(predict_py: Path, deeptmhmm_dir: Path, seq_file: Path, out_dir: Path):
    sample_name = seq_file.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / sample_name

    if dest.exists():
        raise RuntimeError(f'Output directory already exists: {dest}')

    debug_package_layout(deeptmhmm_dir, predict_py)

    cmd = [
        sys.executable,
        str(predict_py),
        '--fasta',
        str(seq_file),
        '--output-dir',
        str(dest),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(deeptmhmm_dir.resolve()),
    )

    if result.returncode != 0:
        raise RuntimeError(
            f'DeepTMHMM failed for {sample_name}\n'
            f'Command: {" ".join(cmd)}\n'
            f'CWD: {deeptmhmm_dir.resolve()}\n'
            f'STDOUT:\n{result.stdout}\n'
            f'STDERR:\n{result.stderr}'
        )

    plot_png = dest / 'plot.png'
    renamed_png = dest / f'{sample_name}.png'
    if plot_png.exists() and not renamed_png.exists():
        plot_png.rename(renamed_png)

    return {
        'sequence_id': sample_name,
        'result_dir': str(dest),
        'png_exists': renamed_png.exists(),
        'success': True,
        'error': '',
    }


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('fasta', type=Path)
    p.add_argument('-o', '--outdir', type=Path, required=True)
    p.add_argument('--deeptmhmm-dir', type=Path, required=True)
    p.add_argument('--manifest', type=Path, default=None)
    p.add_argument('--keep-going', action='store_true')
    args = p.parse_args(argv)

    if not args.fasta.is_file():
        p.error(f'Input FASTA not found: {args.fasta}')
    if not args.deeptmhmm_dir.is_dir():
        p.error(f'DeepTMHMM directory not found: {args.deeptmhmm_dir}')

    predict_py = find_predict_script(args.deeptmhmm_dir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    sequences = [(safe_name(h), s) for h, s in iter_fasta(args.fasta)]
    if not sequences:
        p.error('Input FASTA appears to be empty.')

    rows = []
    tic = time.time()
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        for i, (header, seq) in enumerate(sequences, 1):
            tmp_fasta = tmpdir / f'{header}.fasta'
            with tmp_fasta.open('w') as fh:
                fh.write(f'>{header}\n{seq}\n')
            print(f'[{i}/{len(sequences)}] Running DeepTMHMM on {header}', flush=True)
            try:
                row = run_deeptmhmm_local(predict_py, args.deeptmhmm_dir, tmp_fasta, args.outdir)
            except Exception as e:
                row = {
                    'sequence_id': header,
                    'result_dir': '',
                    'png_exists': False,
                    'success': False,
                    'error': str(e),
                }
                print(f'[ERROR] DeepTMHMM failed for sequence: {header}', file=sys.stderr, flush=True)
                print(str(e), file=sys.stderr, flush=True)
                rows.append(row)
                if not args.keep_going:
                    break
            else:
                rows.append(row)

    manifest = args.manifest or (args.outdir / 'deeptmhmm_manifest.tsv')
    with manifest.open('w', newline='') as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=['sequence_id', 'result_dir', 'png_exists', 'success', 'error'],
            delimiter='\t',
        )
        writer.writeheader()
        writer.writerows(rows)

    failed = [r for r in rows if not r['success']]
    print(f'Processed {len(rows)} sequences in {time.time() - tic:.1f}s', flush=True)
    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
