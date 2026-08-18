#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio.Seq import Seq

REGISTRY_COLUMNS = [
    "internal_id", "source_id", "source_header", "safe_id", "source_type",
    "sequence_sha256", "has_internal_stop", "contig", "start", "end",
    "strand", "gene_id", "transcript_id", "final_output_id", "annotation_track",
    "cds_signature", "duplicate_member_count", "duplicate_members",
]
TRANSLATION_REPORT_COLUMNS = [
    "transcript_id", "status", "reason", "raw_cds_length", "translated_cds_length",
    "trailing_nt_ignored", "internal_stops_replaced", "annotation_track",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize IRnotator inputs into internal-ID proteins and registry.")
    parser.add_argument("--proteins-faa")
    parser.add_argument("--genome-fasta")
    parser.add_argument("--annot-gff")
    parser.add_argument("--out-faa", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--translation-report", default="translation_report.tsv")
    parser.add_argument("--gff-protein-attribute", default="protein_id")
    parser.add_argument("--translation-table", type=int, default=1)
    parser.add_argument("--allow-internal-stops", action="store_true")
    parser.add_argument("--id-prefix", default="IRN")
    parser.add_argument("--annotation-track", default="default")
    parser.add_argument("--start-index", type=int, default=1)
    return parser.parse_args()


def select_mode(args: argparse.Namespace) -> str:
    if args.proteins_faa and not args.genome_fasta and not args.annot_gff:
        return "protein_fasta"
    if not args.proteins_faa and args.genome_fasta and args.annot_gff:
        return "genome_gff"
    if args.proteins_faa and not args.genome_fasta and args.annot_gff:
        return "protein_gff"
    fail("provide proteins FASTA; genome FASTA plus GFF; or proteins FASTA plus GFF")


def read_fasta(path: Path):
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(chunks).replace(" ", "").replace("\t", "").upper()
                    if not sequence:
                        fail(f"empty FASTA sequence for {header!r} in {path}")
                    yield header, sequence
                header = line[1:]
                if not header:
                    fail(f"empty FASTA header at {path}:{line_number}")
                chunks = []
            elif line.strip():
                if header is None:
                    fail(f"FASTA sequence before first header at {path}:{line_number}")
                chunks.append(line.strip())
    if header is None:
        fail(f"no FASTA records found in {path}")
    sequence = "".join(chunks).replace(" ", "").replace("\t", "").upper()
    if not sequence:
        fail(f"empty FASTA sequence for {header!r} in {path}")
    yield header, sequence


def header_id(header: str) -> str:
    value = header.split(maxsplit=1)[0]
    if not value:
        fail(f"cannot derive source_id from FASTA header {header!r}")
    return value


def parse_attributes(raw: str) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"invalid GFF3 attribute without '=': {item!r}")
        key, value = item.split("=", 1)
        if not key:
            raise ValueError("empty GFF3 attribute key")
        result[unquote(key)] = [unquote(part) for part in value.split(",")]
    return result


def parse_gff3(path: Path):
    cds_by_tx: dict[str, list[dict]] = defaultdict(list)
    tx_gene: dict[str, str] = {}
    invalid: dict[str, str] = {}
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue
            contig, _, kind, start, end, _, strand, phase, raw_attrs = fields
            if kind not in {"CDS", "mRNA", "transcript"}:
                continue
            try:
                attributes = parse_attributes(raw_attrs)
            except ValueError as exc:
                if kind == "CDS":
                    invalid.setdefault(f"line:{line_number}", str(exc))
                continue
            feature_id = attributes.get("ID", [None])[0]
            parents = attributes.get("Parent", [])
            if kind in {"mRNA", "transcript"}:
                if feature_id:
                    tx_gene[feature_id] = parents[0] if parents else ""
                continue
            transcript_id = parents[0] if len(parents) == 1 and parents[0] else f"line:{line_number}"
            try:
                start_i, end_i = int(start), int(end)
                if start_i < 1 or end_i < start_i or strand not in {"+", "-"} or phase not in {"0", "1", "2"}:
                    raise ValueError("invalid CDS coordinates, strand, or phase")
            except ValueError as exc:
                invalid.setdefault(transcript_id, str(exc))
                continue
            if len(parents) != 1 or not parents[0]:
                invalid.setdefault(transcript_id, "CDS must have exactly one non-empty Parent attribute")
                continue
            cds_by_tx[transcript_id].append({"contig": contig, "start": start_i, "end": end_i, "strand": strand, "phase": int(phase), "attributes": attributes})
    return cds_by_tx, tx_gene, invalid


def read_genome(path: Path) -> dict[str, str]:
    genome: dict[str, str] = {}
    for header, sequence in read_fasta(path):
        key = header_id(header)
        if key in genome:
            fail(f"duplicate genome contig ID {key!r}")
        genome[key] = sequence
    return genome


def genomic_metadata(transcript_id: str, cds: list[dict], tx_gene: dict[str, str]) -> dict[str, str]:
    contigs = {row["contig"] for row in cds}
    strands = {row["strand"] for row in cds}
    if len(contigs) != 1 or len(strands) != 1:
        raise ValueError("CDS spans multiple contigs or strands")
    ordered = sorted(cds, key=lambda row: (row["start"], row["end"], row["phase"]))
    signature = f"{next(iter(contigs))}|{next(iter(strands))}|" + ",".join(f"{row['start']}-{row['end']}-{row['phase']}" for row in ordered)
    return {"contig": next(iter(contigs)), "start": str(min(row["start"] for row in cds)), "end": str(max(row["end"] for row in cds)), "strand": next(iter(strands)), "gene_id": tx_gene.get(transcript_id, ""), "transcript_id": transcript_id, "cds_signature": signature}


def translate(cds: list[dict], genome: dict[str, str], table: int) -> tuple[str, dict]:
    contig, strand = cds[0]["contig"], cds[0]["strand"]
    if contig not in genome:
        raise ValueError(f"references missing genome contig {contig!r}")
    ordered = sorted(cds, key=lambda row: row["start"], reverse=(strand == "-"))
    fragments = []
    for row in ordered:
        fragment = genome[contig][row["start"] - 1:row["end"]]
        if len(fragment) != row["end"] - row["start"] + 1:
            raise ValueError("CDS coordinates exceed contig bounds")
        fragments.append(str(Seq(fragment).reverse_complement()) if strand == "-" else fragment)
    coding = "".join(fragments)
    trailing = len(coding) % 3
    coding = coding[:-trailing] if trailing else coding
    if not coding:
        raise ValueError("empty CDS after removing incomplete trailing codon")
    protein = str(Seq(coding).translate(table=table, cds=False))
    if protein.endswith("*"):
        protein = protein[:-1]
    stops = protein.count("*")
    protein = protein.replace("*", "X")
    if not protein:
        raise ValueError("empty translated protein")
    return protein, {"raw_cds_length": len("".join(fragments)), "translated_cds_length": len(coding), "trailing_nt_ignored": trailing, "internal_stops_replaced": stops}


def make_row(internal_id: str, source_id: str, source_header: str, source_type: str, sequence: str, track: str, genomic: dict[str, str] | None = None, has_stops: bool = False) -> dict[str, str]:
    genomic = genomic or {}
    return {"internal_id": internal_id, "source_id": source_id, "source_header": source_header, "safe_id": internal_id, "source_type": source_type, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "has_internal_stop": str(has_stops).lower(), "contig": genomic.get("contig", ""), "start": genomic.get("start", ""), "end": genomic.get("end", ""), "strand": genomic.get("strand", ""), "gene_id": genomic.get("gene_id", ""), "transcript_id": genomic.get("transcript_id", ""), "final_output_id": "", "annotation_track": track, "cds_signature": genomic.get("cds_signature", ""), "duplicate_member_count": "1", "duplicate_members": ""}


def main() -> None:
    args = parse_args()
    if args.start_index < 1:
        fail("--start-index must be >= 1")
    mode = select_mode(args)
    records: list[tuple[dict[str, str], str]] = []
    reports: list[dict[str, str]] = []
    if mode == "protein_fasta":
        for index, (header, sequence) in enumerate(read_fasta(Path(args.proteins_faa)), args.start_index):
            records.append((make_row(f"{args.id_prefix}_{index:06d}", header_id(header), header, mode, sequence, args.annotation_track), sequence))
    else:
        cds_by_tx, tx_gene, invalid = parse_gff3(Path(args.annot_gff))
        if mode == "genome_gff":
            genome = read_genome(Path(args.genome_fasta))
            index = args.start_index - 1
            for transcript_id, reason in sorted(invalid.items()):
                reports.append({"transcript_id": transcript_id, "status": "skipped", "reason": reason, "raw_cds_length": "", "translated_cds_length": "", "trailing_nt_ignored": "", "internal_stops_replaced": "", "annotation_track": args.annotation_track})
            for transcript_id in sorted(cds_by_tx):
                if transcript_id in invalid:
                    continue
                try:
                    genomic = genomic_metadata(transcript_id, cds_by_tx[transcript_id], tx_gene)
                    protein, qc = translate(cds_by_tx[transcript_id], genome, args.translation_table)
                except ValueError as exc:
                    print(f"WARNING: skipping {transcript_id}: {exc}", file=sys.stderr)
                    continue
                index += 1
                records.append((make_row(f"{args.id_prefix}_{index:06d}", transcript_id, transcript_id, mode, protein, args.annotation_track, genomic, qc["internal_stops_replaced"] > 0), protein))
                if qc["trailing_nt_ignored"] or qc["internal_stops_replaced"]:
                    reports.append({"transcript_id": transcript_id, "status": "translated_with_warning", "reason": "trailing incomplete codon ignored and/or internal stops replaced with X", "annotation_track": args.annotation_track, **{key: str(value) for key, value in qc.items()}})
        else:
            attributes: dict[str, set[str]] = defaultdict(set)
            for transcript_id, cds in cds_by_tx.items():
                values = {value for row in cds for value in row["attributes"].get(args.gff_protein_attribute, [])}
                if len(values) == 1:
                    attributes[next(iter(values))].add(transcript_id)
            for index, (header, sequence) in enumerate(read_fasta(Path(args.proteins_faa)), args.start_index):
                source_id = header_id(header)
                matches = attributes.get(source_id, set())
                if len(matches) != 1:
                    fail(f"FASTA source_id {source_id!r} has {len(matches)} GFF matches using {args.gff_protein_attribute!r}")
                transcript_id = next(iter(matches))
                records.append((make_row(f"{args.id_prefix}_{index:06d}", source_id, header, mode, sequence, args.annotation_track, genomic_metadata(transcript_id, cds_by_tx[transcript_id], tx_gene)), sequence))
    with Path(args.out_faa).open("w") as fasta:
        for row, sequence in records:
            fasta.write(f">{row['internal_id']}\n")
            for offset in range(0, len(sequence), 60):
                fasta.write(sequence[offset:offset + 60] + "\n")
    with Path(args.registry).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REGISTRY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(row for row, _ in records)
    with Path(args.translation_report).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRANSLATION_REPORT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows({field: row.get(field, "") for field in TRANSLATION_REPORT_COLUMNS} for row in reports)


if __name__ == "__main__":
    main()
