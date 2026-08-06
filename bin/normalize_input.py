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
    "internal_id",
    "source_id",
    "source_header",
    "safe_id",
    "source_type",
    "sequence_sha256",
    "has_internal_stop",
    "contig",
    "start",
    "end",
    "strand",
    "gene_id",
    "transcript_id",
    "final_output_id",
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
    parser.add_argument("--gff-protein-attribute", default="protein_id")
    parser.add_argument("--translation-table", type=int, default=1)
    parser.add_argument("--allow-internal-stops", action="store_true")
    parser.add_argument("--id-prefix", default="IRN")
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
        "provide exactly one valid input combination: "
        "(1) --proteins-faa; "
        "(2) --genome-fasta plus --annot-gff; or "
        "(3) --proteins-faa plus --annot-gff. "
        "Providing all three inputs is not allowed."
    )


def read_fasta(path: Path):
    header = None
    sequence_parts = []

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
    attributes = {}

    for item in raw_attributes.split(";"):
        if not item:
            continue
        if "=" not in item:
            fail(f"invalid GFF3 attribute without '=': {item!r}")

        key, value = item.split("=", 1)
        if not key:
            fail(f"empty GFF3 attribute key in {raw_attributes!r}")

        attributes[unquote(key)] = [unquote(x) for x in value.split(",")]

    return attributes


def parse_gff3(path: Path):
    saw_gff3_header = False
    cds_by_transcript = defaultdict(list)
    transcript_gene = {}

    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")

            if not line:
                continue
            if line == "##gff-version 3":
                saw_gff3_header = True
                continue
            if line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                fail(f"expected 9 GFF3 columns at {path}:{line_number}")

            contig, _, feature_type, start, end, _, strand, phase, raw_attributes = fields

            try:
                start_int = int(start)
                end_int = int(end)
            except ValueError:
                fail(f"invalid GFF3 coordinates at {path}:{line_number}")

            if start_int < 1 or end_int < start_int:
                fail(f"invalid GFF3 interval at {path}:{line_number}")
            if strand not in ("+", "-"):
                fail(f"invalid GFF3 strand at {path}:{line_number}: {strand!r}")

            attributes = parse_attributes(raw_attributes)
            feature_id = attributes.get("ID", [None])[0]
            parents = attributes.get("Parent", [])

            if feature_type in {"mRNA", "transcript"} and feature_id:
                if len(parents) > 1:
                    fail(f"multiple parents for transcript {feature_id!r}")
                transcript_gene[feature_id] = parents[0] if parents else ""

            if feature_type != "CDS":
                continue

            if len(parents) != 1:
                fail(
                    f"CDS at {path}:{line_number} must have exactly one Parent attribute"
                )
            if phase not in ("0", "1", "2"):
                fail(
                    f"CDS at {path}:{line_number} must have phase 0, 1, or 2; "
                    f"got {phase!r}"
                )

            transcript_id = parents[0]
            cds_by_transcript[transcript_id].append(
                {
                    "contig": contig,
                    "start": start_int,
                    "end": end_int,
                    "strand": strand,
                    "phase": int(phase),
                    "attributes": attributes,
                }
            )

    if not saw_gff3_header:
        fail(f"{path} is not strict GFF3: missing '##gff-version 3'")
    if not cds_by_transcript:
        fail(f"no CDS features found in {path}")

    return cds_by_transcript, transcript_gene


def read_genome(path: Path) -> dict[str, str]:
    genome = {}

    for header, sequence in read_fasta(path):
        contig = header_id(header)
        if contig in genome:
            fail(f"duplicate genome contig ID {contig!r} in {path}")
        genome[contig] = sequence

    return genome


def transcript_metadata(transcript_id: str, cds_rows: list[dict], transcript_gene: dict):
    contigs = {row["contig"] for row in cds_rows}
    strands = {row["strand"] for row in cds_rows}

    if len(contigs) != 1:
        fail(f"CDS for transcript {transcript_id!r} spans multiple contigs")
    if len(strands) != 1:
        fail(f"CDS for transcript {transcript_id!r} has inconsistent strands")

    return {
        "contig": next(iter(contigs)),
        "start": str(min(row["start"] for row in cds_rows)),
        "end": str(max(row["end"] for row in cds_rows)),
        "strand": next(iter(strands)),
        "gene_id": transcript_gene.get(transcript_id, ""),
        "transcript_id": transcript_id,
    }


def translate_cds(
    transcript_id: str,
    cds_rows: list[dict],
    genome: dict[str, str],
    translation_table: int,
    allow_internal_stops: bool,
) -> tuple[str, bool]:
    metadata = transcript_metadata(transcript_id, cds_rows, {})
    contig = metadata["contig"]
    strand = metadata["strand"]

    if contig not in genome:
        fail(f"transcript {transcript_id!r} references missing genome contig {contig!r}")

    ordered_rows = sorted(
        cds_rows,
        key=lambda row: row["start"],
        reverse=(strand == "-"),
    )

    parts = []
    for row in ordered_rows:
        fragment = genome[contig][row["start"] - 1 : row["end"]]

        if len(fragment) != row["end"] - row["start"] + 1:
            fail(f"CDS coordinates exceed contig bounds for {transcript_id!r}")

        coding_fragment = str(Seq(fragment).reverse_complement()) if strand == "-" else fragment
        parts.append(coding_fragment[row["phase"] :])

    cds_sequence = "".join(parts)

    if len(cds_sequence) % 3 != 0:
        fail(
            f"translated CDS length is not divisible by 3 for transcript "
            f"{transcript_id!r}"
        )

    try:
        protein = str(Seq(cds_sequence).translate(table=translation_table))
    except Exception as exc:
        fail(f"translation failed for transcript {transcript_id!r}: {exc}")

    if protein.endswith("*"):
        protein = protein[:-1]

    has_internal_stop = "*" in protein
    if has_internal_stop and not allow_internal_stops:
        fail(f"internal stop codon in transcript {transcript_id!r}")

    if has_internal_stop:
        protein = protein.replace("*", "X")

    if not protein:
        fail(f"empty translated protein for transcript {transcript_id!r}")

    return protein, has_internal_stop


def registry_row(
    internal_id: str,
    source_id: str,
    source_header: str,
    source_type: str,
    sequence: str,
    has_internal_stop: bool,
    genomic: dict | None = None,
) -> dict:
    genomic = genomic or {}

    return {
        "internal_id": internal_id,
        "source_id": source_id,
        "source_header": source_header,
        "safe_id": internal_id,
        "source_type": source_type,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "has_internal_stop": str(has_internal_stop).lower(),
        "contig": genomic.get("contig", ""),
        "start": genomic.get("start", ""),
        "end": genomic.get("end", ""),
        "strand": genomic.get("strand", ""),
        "gene_id": genomic.get("gene_id", ""),
        "transcript_id": genomic.get("transcript_id", ""),
        "final_output_id": "",
    }


def write_outputs(records: list[tuple[dict, str]], out_faa: Path, registry: Path) -> None:
    with out_faa.open("w") as fasta_handle:
        for row, sequence in records:
            fasta_handle.write(f">{row['internal_id']}\n")
            for start in range(0, len(sequence), 60):
                fasta_handle.write(f"{sequence[start:start + 60]}\n")

    with registry.open("w", newline="") as registry_handle:
        writer = csv.DictWriter(
            registry_handle,
            fieldnames=REGISTRY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(row for row, _ in records)


def main() -> None:
    args = parse_args()
    mode = select_mode(args)

    out_faa = Path(args.out_faa)
    registry = Path(args.registry)
    records = []

    if mode == "protein_fasta":
        for index, (header, sequence) in enumerate(read_fasta(Path(args.proteins_faa)), start=1):
            internal_id = f"{args.id_prefix}_{index:06d}"
            records.append(
                (
                    registry_row(
                        internal_id=internal_id,
                        source_id=header_id(header),
                        source_header=header,
                        source_type="protein_fasta",
                        sequence=sequence,
                        has_internal_stop=False,
                    ),
                    sequence,
                )
            )

    else:
        cds_by_transcript, transcript_gene = parse_gff3(Path(args.annot_gff))

        if mode == "genome_gff":
            genome = read_genome(Path(args.genome_fasta))

            for index, transcript_id in enumerate(sorted(cds_by_transcript), start=1):
                cds_rows = cds_by_transcript[transcript_id]
                protein, has_internal_stop = translate_cds(
                    transcript_id,
                    cds_rows,
                    genome,
                    args.translation_table,
                    args.allow_internal_stops,
                )
                internal_id = f"{args.id_prefix}_{index:06d}"
                genomic = transcript_metadata(
                    transcript_id,
                    cds_rows,
                    transcript_gene,
                )
                records.append(
                    (
                        registry_row(
                            internal_id=internal_id,
                            source_id=transcript_id,
                            source_header=transcript_id,
                            source_type="genome_gff",
                            sequence=protein,
                            has_internal_stop=has_internal_stop,
                            genomic=genomic,
                        ),
                        protein,
                    )
                )

        else:
            attribute_to_transcripts = defaultdict(set)

            for transcript_id, cds_rows in cds_by_transcript.items():
                values = set()
                for cds_row in cds_rows:
                    values.update(
                        cds_row["attributes"].get(args.gff_protein_attribute, [])
                    )

                if len(values) != 1:
                    fail(
                        f"transcript {transcript_id!r} must have exactly one "
                        f"{args.gff_protein_attribute!r} value across CDS features"
                    )

                attribute_to_transcripts[next(iter(values))].add(transcript_id)

            for index, (header, sequence) in enumerate(
                read_fasta(Path(args.proteins_faa)),
                start=1,
            ):
                source_id = header_id(header)
                matches = attribute_to_transcripts.get(source_id, set())

                if len(matches) != 1:
                    fail(
                        f"FASTA source_id {source_id!r} has {len(matches)} GFF3 matches "
                        f"using attribute {args.gff_protein_attribute!r}"
                    )

                transcript_id = next(iter(matches))
                internal_id = f"{args.id_prefix}_{index:06d}"
                genomic = transcript_metadata(
                    transcript_id,
                    cds_by_transcript[transcript_id],
                    transcript_gene,
                )
                records.append(
                    (
                        registry_row(
                            internal_id=internal_id,
                            source_id=source_id,
                            source_header=header,
                            source_type="protein_gff",
                            sequence=sequence,
                            has_internal_stop=False,
                            genomic=genomic,
                        ),
                        sequence,
                    )
                )

    write_outputs(records, out_faa, registry)


if __name__ == "__main__":
    main()
