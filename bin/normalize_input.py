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
]

TRANSLATION_REPORT_COLUMNS = [
    "transcript_id", "status", "reason", "raw_cds_length",
    "translated_cds_length", "trailing_nt_ignored", "internal_stops_replaced",
    "annotation_track",
]


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize IRnotator inputs into internal-ID proteins and registry."
    )
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
    proteins = bool(args.proteins_faa)
    genome = bool(args.genome_fasta)
    gff = bool(args.annot_gff)
    if proteins and not genome and not gff:
        return "protein_fasta"
    if not proteins and genome and gff:
        return "genome_gff"
    if proteins and not genome and gff:
        return "protein_gff"
    fail(
        "provide exactly one valid input combination: (1) --proteins-faa; "
        "(2) --genome-fasta plus --annot-gff; or (3) --proteins-faa plus --annot-gff"
    )


def read_fasta(path: Path):
    header = None
    sequence_parts: list[str] = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_parts).replace(" ", "").replace("\t", "")
                    if not sequence:
                        fail(f"empty FASTA sequence for header {header!r} in {path}")
                    yield header, sequence.upper()
                header = line[1:]
                if not header:
                    fail(f"empty FASTA header at {path}:{line_number}")
                sequence_parts = []
            elif line.strip():
                if header is None:
                    fail(f"FASTA sequence before first header at {path}:{line_number}")
                sequence_parts.append(line.strip())
    if header is None:
        fail(f"no FASTA records found in {path}")
    sequence = "".join(sequence_parts).replace(" ", "").replace("\t", "")
    if not sequence:
        fail(f"empty FASTA sequence for header {header!r} in {path}")
    yield header, sequence.upper()


def header_id(header: str) -> str:
    source_id = header.split(maxsplit=1)[0]
    if not source_id:
        fail(f"cannot derive source_id from FASTA header {header!r}")
    return source_id


def parse_attributes(raw_attributes: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = {}
    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"invalid GFF3 attribute without '=': {item!r}")
        key, value = item.split("=", 1)
        if not key:
            raise ValueError("empty GFF3 attribute key")
        attributes[unquote(key)] = [unquote(x) for x in value.split(",")]
    return attributes


def parse_gff3(path: Path):
    """Parse only records relevant to protein translation.

    Non-CDS records are ignored except mRNA/transcript records, which are used
    opportunistically to retain gene IDs.  Malformed CDS records are associated
    with their parent transcript where possible and reported as skipped later.
    """
    cds_by_transcript: dict[str, list[dict]] = defaultdict(list)
    transcript_gene: dict[str, str] = {}
    invalid_transcripts: dict[str, str] = {}

    def mark_invalid(transcript_id: str, reason: str) -> None:
        invalid_transcripts.setdefault(transcript_id, reason)

    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                # Irrelevant malformed rows must not stop normalization. A CDS
                # cannot be identified reliably without nine columns.
                continue

            contig, _, feature_type, start, end, _, strand, phase, raw_attributes = fields
            if feature_type not in {"CDS", "mRNA", "transcript"}:
                continue

            try:
                attributes = parse_attributes(raw_attributes)
            except ValueError as exc:
                if feature_type == "CDS":
                    mark_invalid(f"line:{line_number}", str(exc))
                continue

            feature_id = attributes.get("ID", [None])[0]
            parents = attributes.get("Parent", [])
            if feature_type in {"mRNA", "transcript"}:
                if feature_id:
                    transcript_gene[feature_id] = parents[0] if parents else ""
                continue

            transcript_id = parents[0] if len(parents) == 1 and parents[0] else f"line:{line_number}"
            if len(parents) != 1 or not parents[0]:
                mark_invalid(transcript_id, "CDS must have exactly one non-empty Parent attribute")
                continue
            try:
                start_int = int(start)
                end_int = int(end)
            except ValueError:
                mark_invalid(transcript_id, "CDS has non-numeric coordinates")
                continue
            if start_int < 1 or end_int < start_int:
                mark_invalid(transcript_id, "CDS has an invalid genomic interval")
                continue
            if strand not in {"+", "-"}:
                mark_invalid(transcript_id, f"CDS has invalid strand {strand!r}")
                continue
            if phase not in {"0", "1", "2"}:
                mark_invalid(transcript_id, f"CDS has invalid phase {phase!r}")
                continue
            cds_by_transcript[transcript_id].append({
                "contig": contig,
                "start": start_int,
                "end": end_int,
                "strand": strand,
                "phase": int(phase),
                "attributes": attributes,
            })

    return cds_by_transcript, transcript_gene, invalid_transcripts


def read_genome(path: Path) -> dict[str, str]:
    genome: dict[str, str] = {}
    for header, sequence in read_fasta(path):
        contig = header_id(header)
        if contig in genome:
            fail(f"duplicate genome contig ID {contig!r} in {path}")
        genome[contig] = sequence
    return genome


def transcript_metadata(transcript_id: str, cds_rows: list[dict], transcript_gene: dict) -> dict:
    contigs = {row["contig"] for row in cds_rows}
    strands = {row["strand"] for row in cds_rows}
    if len(contigs) != 1:
        raise ValueError(f"CDS spans multiple contigs: {sorted(contigs)}")
    if len(strands) != 1:
        raise ValueError(f"CDS has inconsistent strands: {sorted(strands)}")
    return {
        "contig": next(iter(contigs)),
        "start": str(min(row["start"] for row in cds_rows)),
        "end": str(max(row["end"] for row in cds_rows)),
        "strand": next(iter(strands)),
        "gene_id": transcript_gene.get(transcript_id, ""),
        "transcript_id": transcript_id,
    }


def translate_cds(cds_rows: list[dict], genome: dict[str, str], translation_table: int) -> tuple[str, dict]:
    contig = cds_rows[0]["contig"]
    strand = cds_rows[0]["strand"]
    if contig not in genome:
        raise ValueError(f"references missing genome contig {contig!r}")
    ordered_rows = sorted(cds_rows, key=lambda row: row["start"], reverse=(strand == "-"))
    parts: list[str] = []
    for row in ordered_rows:
        fragment = genome[contig][row["start"] - 1:row["end"]]
        if len(fragment) != row["end"] - row["start"] + 1:
            raise ValueError("CDS coordinates exceed contig bounds")
        parts.append(str(Seq(fragment).reverse_complement()) if strand == "-" else fragment)
    cds_sequence = "".join(parts)
    raw_length = len(cds_sequence)
    trailing_nt_ignored = raw_length % 3
    translated_cds = cds_sequence[:raw_length - trailing_nt_ignored] if trailing_nt_ignored else cds_sequence
    if not translated_cds:
        raise ValueError("empty CDS after removing incomplete trailing codon")
    protein = str(Seq(translated_cds).translate(table=translation_table, cds=False))
    if protein.endswith("*"):
        protein = protein[:-1]
    internal_stops_replaced = protein.count("*")
    protein = protein.replace("*", "X")
    if not protein:
        raise ValueError("empty translated protein")
    return protein, {
        "raw_cds_length": raw_length,
        "translated_cds_length": len(translated_cds),
        "trailing_nt_ignored": trailing_nt_ignored,
        "internal_stops_replaced": internal_stops_replaced,
    }


def registry_row(internal_id: str, source_id: str, source_header: str, source_type: str,
                 sequence: str, has_internal_stop: bool, annotation_track: str,
                 genomic: dict | None = None) -> dict:
    genomic = genomic or {}
    return {
        "internal_id": internal_id, "source_id": source_id, "source_header": source_header,
        "safe_id": internal_id, "source_type": source_type,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "has_internal_stop": str(has_internal_stop).lower(),
        "contig": genomic.get("contig", ""), "start": genomic.get("start", ""),
        "end": genomic.get("end", ""), "strand": genomic.get("strand", ""),
        "gene_id": genomic.get("gene_id", ""), "transcript_id": genomic.get("transcript_id", ""),
        "final_output_id": "", "annotation_track": annotation_track,
    }


def write_outputs(records: list[tuple[dict, str]], out_faa: Path, registry: Path) -> None:
    with out_faa.open("w") as fasta_handle:
        for row, sequence in records:
            fasta_handle.write(f">{row['internal_id']}\n")
            for start in range(0, len(sequence), 60):
                fasta_handle.write(f"{sequence[start:start + 60]}\n")
    with registry.open("w", newline="") as registry_handle:
        writer = csv.DictWriter(registry_handle, fieldnames=REGISTRY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(row for row, _ in records)


def write_translation_report(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRANSLATION_REPORT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def report_row(transcript_id: str, status: str, reason: str, annotation_track: str, **qc) -> dict:
    return {
        "transcript_id": transcript_id, "status": status, "reason": reason,
        "raw_cds_length": qc.get("raw_cds_length", ""),
        "translated_cds_length": qc.get("translated_cds_length", ""),
        "trailing_nt_ignored": qc.get("trailing_nt_ignored", ""),
        "internal_stops_replaced": qc.get("internal_stops_replaced", ""),
        "annotation_track": annotation_track,
    }


def main() -> None:
    args = parse_args()
    if args.start_index < 1:
        fail("--start-index must be >= 1")
    mode = select_mode(args)
    records: list[tuple[dict, str]] = []
    report_rows: list[dict] = []

    if mode == "protein_fasta":
        for index, (header, sequence) in enumerate(read_fasta(Path(args.proteins_faa)), start=args.start_index):
            internal_id = f"{args.id_prefix}_{index:06d}"
            records.append((registry_row(internal_id, header_id(header), header, "protein_fasta", sequence, False, args.annotation_track), sequence))
    else:
        cds_by_transcript, transcript_gene, invalid_transcripts = parse_gff3(Path(args.annot_gff))
        if mode == "genome_gff":
            genome = read_genome(Path(args.genome_fasta))
            output_index = args.start_index - 1
            for transcript_id, reason in sorted(invalid_transcripts.items()):
                report_rows.append(report_row(transcript_id, "skipped", reason, args.annotation_track))
            for transcript_id in sorted(cds_by_transcript):
                if transcript_id in invalid_transcripts:
                    continue
                try:
                    cds_rows = cds_by_transcript[transcript_id]
                    genomic = transcript_metadata(transcript_id, cds_rows, transcript_gene)
                    protein, qc = translate_cds(cds_rows, genome, args.translation_table)
                except Exception as exc:
                    message = str(exc)
                    report_rows.append(report_row(transcript_id, "skipped", message, args.annotation_track))
                    print(f"WARNING: skipping {transcript_id}: {message}", file=sys.stderr)
                    continue
                output_index += 1
                internal_id = f"{args.id_prefix}_{output_index:06d}"
                records.append((registry_row(internal_id, transcript_id, transcript_id, "genome_gff", protein, qc["internal_stops_replaced"] > 0, args.annotation_track, genomic), protein))
                if qc["trailing_nt_ignored"] or qc["internal_stops_replaced"]:
                    report_rows.append(report_row(transcript_id, "translated_with_warning", "trailing incomplete codon ignored and/or internal stops replaced with X", args.annotation_track, **qc))
        else:
            attribute_to_transcripts: dict[str, set[str]] = defaultdict(set)
            for transcript_id, cds_rows in cds_by_transcript.items():
                values = {value for row in cds_rows for value in row["attributes"].get(args.gff_protein_attribute, [])}
                if len(values) == 1:
                    attribute_to_transcripts[next(iter(values))].add(transcript_id)
            for index, (header, sequence) in enumerate(read_fasta(Path(args.proteins_faa)), start=args.start_index):
                source_id = header_id(header)
                matches = attribute_to_transcripts.get(source_id, set())
                if len(matches) != 1:
                    fail(f"FASTA source_id {source_id!r} has {len(matches)} GFF matches using {args.gff_protein_attribute!r}")
                transcript_id = next(iter(matches))
                genomic = transcript_metadata(transcript_id, cds_by_transcript[transcript_id], transcript_gene)
                internal_id = f"{args.id_prefix}_{index:06d}"
                records.append((registry_row(internal_id, source_id, header, "protein_gff", sequence, False, args.annotation_track, genomic), sequence))

    write_outputs(records, Path(args.out_faa), Path(args.registry))
    write_translation_report(report_rows, Path(args.translation_report))


if __name__ == "__main__":
    main()
