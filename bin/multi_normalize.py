#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def fasta_records(path: Path):
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
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
    if header is None or not chunks:
        raise ValueError(f"no complete FASTA records in {path}")
    yield header, "".join(chunks)


def gff_seqids(path: Path) -> set[str]:
    result = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) == 9:
                result.add(fields[0])
    return result


def track_name(path: Path) -> str:
    return re.sub(r"\.(gff3?|gtf)$", "", path.name, flags=re.IGNORECASE)


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), reader.fieldnames or []


def write_fasta(records: list[tuple[dict[str, str], str]], path: Path) -> None:
    with path.open("w") as handle:
        for row, sequence in records:
            handle.write(f">{row['internal_id']}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset:offset + 60] + "\n")


def provenance_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (row["annotation_track"], row["source_id"], row["transcript_id"])


def member_label(row: dict[str, str]) -> str:
    return ":".join(provenance_key(row))


def canonicalize(records: list[tuple[dict[str, str], str]]) -> tuple[list[tuple[dict[str, str], str]], list[dict[str, str]]]:
    by_signature: dict[str, list[tuple[dict[str, str], str]]] = defaultdict(list)
    for row, sequence in records:
        signature = row.get("cds_signature", "")
        key = signature if signature else f"unique:{row['internal_id']}"
        by_signature[key].append((row, sequence))

    grouped = []
    for key, members in by_signature.items():
        members.sort(key=lambda record: provenance_key(record[0]))
        signature = members[0][0].get("cds_signature", "")
        if signature:
            hashes = {row["sequence_sha256"] for row, _ in members}
            if len(hashes) != 1:
                detail = "; ".join(member_label(row) for row, _ in members)
                raise ValueError(f"exact CDS duplicate group has non-identical translations: {signature}; members: {detail}")
        grouped.append(members)

    grouped.sort(key=lambda members: provenance_key(members[0][0]))
    canonical_records: list[tuple[dict[str, str], str]] = []
    provenance: list[dict[str, str]] = []
    for index, members in enumerate(grouped, 1):
        canonical, sequence = members[0]
        canonical_id = f"IRN_{index:06d}"
        labels = [member_label(row) for row, _ in members]
        canonical = dict(canonical)
        canonical["internal_id"] = canonical_id
        canonical["safe_id"] = canonical_id
        canonical["duplicate_member_count"] = str(len(members))
        canonical["duplicate_members"] = ",".join(labels)
        canonical_records.append((canonical, sequence))
        for row, _ in members:
            provenance.append({
                "canonical_internal_id": canonical_id,
                "provisional_internal_id": row["internal_id"],
                "annotation_track": row["annotation_track"],
                "source_id": row["source_id"],
                "transcript_id": row["transcript_id"],
                "cds_signature": row.get("cds_signature", ""),
                "sequence_sha256": row["sequence_sha256"],
            })
    return canonical_records, provenance


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize multi-GFF genome annotation and canonicalize exact CDS duplicates.")
    parser.add_argument("--genome-fasta", nargs="+", required=True)
    parser.add_argument("--annot-gff", nargs="+", required=True)
    parser.add_argument("--normalizer", required=True)
    parser.add_argument("--out-faa", required=True)
    parser.add_argument("--out-registry", required=True)
    parser.add_argument("--out-report", required=True)
    parser.add_argument("--out-duplicate-provenance", required=True)
    parser.add_argument("--gff-protein-attribute", default="protein_id")
    parser.add_argument("--translation-table", type=int, default=1)
    parser.add_argument("--allow-internal-stops", action="store_true")
    args = parser.parse_args()

    genomes = [Path(value) for value in args.genome_fasta]
    gffs = [Path(value) for value in args.annot_gff]
    normalizer = Path(args.normalizer)
    for path in [*genomes, *gffs, normalizer]:
        if not path.exists():
            raise ValueError(f"input file does not exist: {path}")
        if path.suffix.lower() == ".gz":
            raise ValueError(f"compressed inputs are not supported: {path}")

    tracks = [track_name(path) for path in gffs]
    if not all(tracks) or len(set(tracks)) != len(tracks):
        raise ValueError("annotation GFF filenames must produce unique non-empty track names")

    with tempfile.TemporaryDirectory() as temp_text:
        temp = Path(temp_text)
        genome_path = temp / "genome.fna"
        seen_contigs = set()
        with genome_path.open("w") as output:
            for source in genomes:
                for seqid, sequence in fasta_records(source):
                    if seqid in seen_contigs:
                        raise ValueError(f"duplicate FASTA sequence ID is forbidden: {seqid!r}")
                    seen_contigs.add(seqid)
                    output.write(f">{seqid}\n{sequence}\n")
        for gff in gffs:
            missing = gff_seqids(gff) - seen_contigs
            if missing:
                raise ValueError(f"{gff}: GFF seqid absent from accumulated FASTA: {sorted(missing)[0]!r}")

        all_records: list[tuple[dict[str, str], str]] = []
        report_rows: list[dict[str, str]] = []
        registry_fields: list[str] | None = None
        provisional_index = 1
        for gff, track in zip(gffs, tracks):
            track_dir = temp / track
            track_dir.mkdir()
            track_faa = track_dir / "proteins.faa"
            track_registry = track_dir / "registry.tsv"
            track_report = track_dir / "report.tsv"
            command = [sys.executable, str(normalizer), "--genome-fasta", str(genome_path), "--annot-gff", str(gff), "--annotation-track", track, "--start-index", str(provisional_index), "--out-faa", str(track_faa), "--registry", str(track_registry), "--translation-report", str(track_report), "--gff-protein-attribute", args.gff_protein_attribute, "--translation-table", str(args.translation_table)]
            if args.allow_internal_stops:
                command.append("--allow-internal-stops")
            subprocess.run(command, check=True)
            rows, fields = read_tsv(track_registry)
            if registry_fields is None:
                registry_fields = fields
            elif fields != registry_fields:
                raise ValueError(f"registry schema differs between annotation tracks: {gff}")
            sequences = dict(fasta_records(track_faa))
            if set(sequences) != {row["internal_id"] for row in rows}:
                raise ValueError(f"FASTA/registry ID mismatch after normalizing {gff}")
            all_records.extend((row, sequences[row["internal_id"]]) for row in rows)
            report_rows.extend(read_tsv(track_report)[0])
            provisional_index += len(rows)

    if not registry_fields:
        raise ValueError("no annotation GFF records were normalized")
    canonical_records, provenance = canonicalize(all_records)
    write_fasta(canonical_records, Path(args.out_faa))
    with Path(args.out_registry).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=registry_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(row for row, _ in canonical_records)
    with Path(args.out_duplicate_provenance).open("w", newline="") as handle:
        fields = ["canonical_internal_id", "provisional_internal_id", "annotation_track", "source_id", "transcript_id", "cds_signature", "sequence_sha256"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(provenance)
    with Path(args.out_report).open("w", newline="") as handle:
        fields = ["transcript_id", "annotation_track", "status", "reason", "raw_cds_length", "translated_cds_length", "trailing_nt_ignored", "internal_stops_replaced"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows({field: row.get(field, "") for field in fields} for row in report_rows)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
