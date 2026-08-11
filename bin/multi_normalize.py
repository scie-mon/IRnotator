#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def fasta_records(path: Path):
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    if not chunks:
                        raise ValueError(f"empty FASTA sequence for {header!r} in {path}")
                    yield header, "".join(chunks)
                header = line[1:].split(maxsplit=1)[0]
                if not header:
                    raise ValueError(f"empty FASTA header at {path}:{line_number}")
                chunks = []
            elif line.strip():
                if header is None:
                    raise ValueError(f"FASTA sequence before header at {path}:{line_number}")
                chunks.append(line.strip())
    if header is None:
        raise ValueError(f"no FASTA records in {path}")
    if not chunks:
        raise ValueError(f"empty FASTA sequence for {header!r} in {path}")
    yield header, "".join(chunks)


def gff_seqids(path: Path) -> set[str]:
    """Collect seqids from structurally GFF-like feature rows only.

    Non-feature metadata and malformed non-CDS rows are deliberately left to
    the single-track normalizer, which can ignore them or report affected CDS
    records without aborting unrelated annotation tracks.
    """
    seqids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue
            seqids.add(fields[0])
    return seqids


def track_name(path: Path) -> str:
    return re.sub(r"\.(gff3?|gtf)$", "", path.name, flags=re.IGNORECASE)


def append_file(source: Path, destination: Path) -> None:
    with source.open("r") as source_handle, destination.open("a") as destination_handle:
        shutil.copyfileobj(source_handle, destination_handle)


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames or []


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize multiple genome/GFF tracks sequentially into global IDs."
    )
    parser.add_argument("--genome-fasta", nargs="+", required=True)
    parser.add_argument("--annot-gff", nargs="+", required=True)
    parser.add_argument("--normalizer", required=True)
    parser.add_argument("--out-faa", required=True)
    parser.add_argument("--out-registry", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--gff-protein-attribute", default="protein_id")
    parser.add_argument("--translation-table", type=int, default=1)
    parser.add_argument("--allow-internal-stops", action="store_true")
    args = parser.parse_args()

    genome_paths = [Path(path) for path in args.genome_fasta]
    gff_paths = [Path(path) for path in args.annot_gff]
    normalizer_path = Path(args.normalizer)

    for path in [*genome_paths, *gff_paths, normalizer_path]:
        if not path.exists():
            raise ValueError(f"input file does not exist: {path}")
        if path.suffix.lower() == ".gz":
            raise ValueError(f"compressed inputs are not supported: {path}")

    genome_ids: set[str] = set()
    track_names: set[str] = set()
    for gff in gff_paths:
        current_track = track_name(gff)
        if not current_track:
            raise ValueError(f"cannot derive annotation track name from {gff}")
        if current_track in track_names:
            raise ValueError(
                "duplicate annotation track name derived from GFF filename: "
                f"{current_track!r}"
            )
        track_names.add(current_track)

    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_directory = Path(temporary_directory)
        accumulated_genome = temporary_directory / "accumulated_genome.fna"
        with accumulated_genome.open("w") as output:
            for fasta in genome_paths:
                for seqid, sequence in fasta_records(fasta):
                    if seqid in genome_ids:
                        raise ValueError(f"duplicate FASTA sequence ID is forbidden: {seqid!r}")
                    genome_ids.add(seqid)
                    output.write(f">{seqid}\n")
                    for start in range(0, len(sequence), 80):
                        output.write(sequence[start:start + 80] + "\n")

        for gff in gff_paths:
            missing = gff_seqids(gff) - genome_ids
            if missing:
                raise ValueError(
                    f"{gff}: GFF seqid absent from accumulated FASTA: {sorted(missing)[0]!r}"
                )

        out_faa = Path(args.out_faa)
        out_registry = Path(args.out_registry)
        out_report = Path(args.out_report)
        out_faa.write_text("")
        report_rows: list[dict[str, str]] = []
        registry_fields: list[str] | None = None
        global_index = 1

        for gff in gff_paths:
            current_track = track_name(gff)
            track_directory = temporary_directory / current_track
            track_directory.mkdir()
            track_faa = track_directory / "normalized_proteins.faa"
            track_registry = track_directory / "sequence_registry.tsv"
            track_report = track_directory / "translation_report.tsv"

            command = [
                sys.executable, str(normalizer_path),
                "--genome-fasta", str(accumulated_genome),
                "--annot-gff", str(gff),
                "--annotation-track", current_track,
                "--start-index", str(global_index),
                "--out-faa", str(track_faa),
                "--registry", str(track_registry),
                "--translation-report", str(track_report),
                "--gff-protein-attribute", args.gff_protein_attribute,
                "--translation-table", str(args.translation_table),
            ]
            if args.allow_internal_stops:
                command.append("--allow-internal-stops")
            subprocess.run(command, check=True)

            append_file(track_faa, out_faa)
            track_rows, track_fields = read_tsv(track_registry)
            if registry_fields is None:
                registry_fields = track_fields
                with out_registry.open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(track_rows)
            else:
                if track_fields != registry_fields:
                    raise ValueError(f"registry schema differs between annotation tracks: {gff}")
                with out_registry.open("a", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t", lineterminator="\n")
                    writer.writerows(track_rows)

            if track_report.exists():
                rows, _ = read_tsv(track_report)
                report_rows.extend(rows)
            global_index += len(track_rows)

        if registry_fields is None:
            raise ValueError("no annotation GFF files were supplied")

        report_fields = [
            "transcript_id", "annotation_track", "status", "reason",
            "raw_cds_length", "translated_cds_length", "trailing_nt_ignored",
            "internal_stops_replaced",
        ]
        with out_report.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=report_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows({field: row.get(field, "") for field in report_fields} for row in report_rows)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
