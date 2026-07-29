#!/usr/bin/env python3

import argparse
import csv
import re
import sys
import tempfile
import time
from pathlib import Path

import biolib


def iter_fasta(path: Path):
    header = None
    seq_lines = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, "".join(seq_lines)


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def run_deeptmhmm(app, seq_file: Path, out_dir: Path, backend: str):
    args = f"--fasta {seq_file}"
    if backend == "local":
        job = app.cli(args=args, machine="local")
    else:
        job = app.cli(args=args)

    sample_name = seq_file.stem
    dest = out_dir / sample_name
    dest.mkdir(parents=True, exist_ok=True)
    job.save_files(str(dest))

    plot_png = dest / "plot.png"
    renamed_png = dest / f"{sample_name}.png"
    if plot_png.exists() and not renamed_png.exists():
        plot_png.rename(renamed_png)

    return {
        "sequence_id": sample_name,
        "result_dir": str(dest),
        "png_exists": renamed_png.exists(),
        "success": True,
        "error": "",
    }


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="deeptmhmm_runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run DeepTMHMM once per sequence from a multi-FASTA."
    )
    p.add_argument("fasta", type=Path, help="Input FASTA file.")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory.")
    p.add_argument("--backend", choices=["cloud", "local"], default="cloud",
                   help="Execution backend for BioLib.")
    p.add_argument("--app", default="DTU/DeepTMHMM:1.0.24",
                   help="BioLib app identifier.")
    p.add_argument("--manifest", type=Path, default=None,
                   help="Optional TSV manifest path.")
    p.add_argument("--keep-going", action="store_true",
                   help="Continue if one sequence fails.")
    args = p.parse_args(argv)

    if not args.fasta.is_file():
        p.error(f"Input FASTA not found: {args.fasta}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    sequences = [(safe_name(h), s) for h, s in iter_fasta(args.fasta)]
    if not sequences:
        p.error("Input FASTA appears to be empty.")

    app = biolib.load(args.app)
    rows = []
    tic = time.time()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        for i, (header, seq) in enumerate(sequences, 1):
            tmp_fasta = tmpdir / f"{header}.fasta"
            with tmp_fasta.open("w") as fh:
                fh.write(f">{header}\n{seq}\n")

            print(f"[{i}/{len(sequences)}] Running DeepTMHMM on {header}", flush=True)

            try:
                row = run_deeptmhmm(app, tmp_fasta, args.outdir, args.backend)
            except Exception as e:
                row = {
                    "sequence_id": header,
                    "result_dir": "",
                    "png_exists": False,
                    "success": False,
                    "error": str(e),
                }
                if not args.keep_going:
                    rows.append(row)
                    break
            rows.append(row)

    manifest = args.manifest or (args.outdir / "deeptmhmm_manifest.tsv")
    with manifest.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sequence_id", "result_dir", "png_exists", "success", "error"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    failed = [r for r in rows if not r["success"]]
    print(f"Processed {len(rows)} sequences in {time.time() - tic:.1f}s", flush=True)

    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
    