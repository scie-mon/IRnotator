#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import math
import sys
from collections import defaultdict
from pathlib import Path

MARGIN_BP = 20
STATUS = {
    "accept": ("#31a354", "#176b37"),
    "reject": ("#de2d26", "#9f1b17"),
    "review": ("#8c8c8c", "#555555"),
}
NO_SIGNAL = "#f28e2b"
VALID_DECISIONS = set(STATUS)


def attrs(text):
    return {
        key: value
        for item in text.rstrip(";").split(";")
        if "=" in item
        for key, value in [item.split("=", 1)]
    }


def gff(path):
    output = []
    for line_number, line in enumerate(path.open(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"{path}:{line_number}: expected 9-column GFF3")
        output.append(
            dict(
                seqid=fields[0],
                kind=fields[2],
                start=int(fields[3]),
                end=int(fields[4]),
                strand=fields[6],
                attrs=attrs(fields[8]),
            )
        )
    return output


def signals(path):
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not rows.fieldnames or not {"protein_id", "has_signal"}.issubset(rows.fieldnames):
            raise ValueError(f"{path}: requires protein_id and has_signal")
        result = {}
        for row in rows:
            value = row["has_signal"].strip().upper()
            if value not in {"TRUE", "FALSE"}:
                raise ValueError(f"{path}: invalid has_signal {value!r}")
            result[row["protein_id"].strip()] = value == "TRUE"
    return result


def decisions(path):
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not rows.fieldnames or not {"protein_id", "final_decision"}.issubset(rows.fieldnames):
            raise ValueError(f"{path}: requires protein_id and final_decision")
        result = {}
        for row in rows:
            protein_id = row["protein_id"].strip()
            decision = row["final_decision"].strip().lower()
            if not protein_id:
                raise ValueError(f"{path}: empty protein_id")
            if decision not in VALID_DECISIONS:
                raise ValueError(f"{path}: invalid final_decision for {protein_id}: {decision!r}")
            if protein_id in result:
                raise ValueError(f"{path}: duplicate protein_id: {protein_id}")
            result[protein_id] = decision
    return result


def pale(colour, fraction=0.68):
    return "#" + "".join(
        f"{round(int(colour[index:index + 2], 16) + (255 - int(colour[index:index + 2], 16)) * fraction):02x}"
        for index in (1, 3, 5)
    )


def poly(points, fill, stroke, width):
    return '<polygon points="{}" fill="{}" stroke="{}" stroke-width="{}"/>'.format(
        " ".join(f"{x:.2f},{y:.2f}" for x, y in points), fill, stroke, width
    )


def block(x1, x2, y, height, strand, fill, stroke, width):
    delta = min(14, max(5, (x2 - x1) / 3))
    points = (
        [(x1, y + height / 2), (x1 + delta, y), (x2, y), (x2, y + height), (x1 + delta, y + height)]
        if strand == "-"
        else [(x1, y), (x2 - delta, y), (x2, y + height / 2), (x2 - delta, y + height), (x1, y + height)]
    )
    return poly(points, fill, stroke, width)


def step(span):
    raw = max(span / 6, 1)
    base = 10 ** math.floor(math.log10(raw))
    return int(next(value * base for value in (1, 2, 5, 10) if raw <= value * base))


def render(gff_path, feature_path, decision_path, out_path, focal_internal):
    records = gff(gff_path)
    signal = signals(feature_path)
    decision_by_id = decisions(decision_path)
    transcripts = {
        record["attrs"].get("ID"): record
        for record in records
        if record["kind"] in {"mRNA", "transcript"} and record["attrs"].get("ID")
    }
    focal = next(
        (
            transcript_id
            for transcript_id, record in transcripts.items()
            if record["attrs"].get("internal_id") == focal_internal
        ),
        None,
    )
    if focal is None:
        raise ValueError(f"focal internal_id absent: {focal_internal}")

    parents = [parent for parent in transcripts[focal]["attrs"].get("Parent", "").split(",") if parent]
    if len(parents) != 1:
        raise ValueError("focal transcript lacks one parent gene")
    gene_id = parents[0]
    genes = [record for record in records if record["kind"] == "gene" and record["attrs"].get("ID") == gene_id]
    if len(genes) != 1:
        raise ValueError(f"expected one gene {gene_id}")
    gene = genes[0]

    transcripts = {
        transcript_id: record
        for transcript_id, record in transcripts.items()
        if gene_id in record["attrs"].get("Parent", "").split(",")
    }
    cds = defaultdict(list)
    for record in records:
        if record["kind"] == "CDS":
            for parent in record["attrs"].get("Parent", "").split(","):
                if parent in transcripts:
                    cds[parent].append(record)
    if set(transcripts) - set(cds):
        raise ValueError("missing CDS records")

    transcript_ids = sorted(
        transcripts,
        key=lambda transcript_id: (
            transcripts[transcript_id]["start"],
            transcripts[transcript_id]["end"],
            transcripts[transcript_id]["attrs"].get("internal_id", ""),
            transcript_id,
        ),
    )
    low, high = gene["start"] - MARGIN_BP, gene["end"] + MARGIN_BP
    span = high - low + 1
    width, left, meta, right = 1200, 110, 190, 35
    draw = width - left - meta - right
    axis, gene_y, row0, row_height, cds_height = 52, 90, 138, 47, 20
    height = row0 + len(transcript_ids) * row_height + 28
    x = lambda base: left + (base - low) * draw / span

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:16px;font-weight:bold}.sub{font-size:11px;fill:#555}.meta{font-size:11px}.tick{font-size:10px;fill:#555}</style>',
        f'<text x="{left}" y="20" class="title">{html.escape(gene_id)}</text>',
        f'<text x="{left}" y="37" class="sub">{html.escape(gene["seqid"])} · {gene["start"]:,}–{gene["end"]:,} · {gene["strand"]} strand</text>',
    ]
    major = step(span)
    minor = major / 5
    tick = math.ceil(low / minor) * minor
    while tick <= high:
        position = x(round(tick))
        is_major = abs(tick / major - round(tick / major)) < 1e-9
        colour = "#cccccc" if is_major else "#e8e8e8"
        svg.append(f'<line x1="{position:.2f}" y1="{axis}" x2="{position:.2f}" y2="{height - 10}" stroke="{colour}"/>')
        if is_major:
            svg.append(f'<text x="{position:.2f}" y="{axis - 5}" text-anchor="middle" class="tick">{int(tick):,}</text>')
        tick += minor

    svg += [
        f'<text x="{left - 8}" y="{gene_y + 10}" text-anchor="end" class="meta">gene</text>',
        f'<text x="{left + draw + 8}" y="{gene_y + 10}" class="meta">{gene["end"] - gene["start"] + 1:,} bp</text>',
        block(x(gene["start"]), x(gene["end"] + 1), gene_y, 12, gene["strand"], "#303030", "#303030", 1),
    ]

    for index, transcript_id in enumerate(transcript_ids):
        transcript = transcripts[transcript_id]
        segments = sorted(cds[transcript_id], key=lambda record: record["start"])
        y = row0 + index * row_height
        internal_id = transcript["attrs"].get("internal_id", transcript_id)
        focal_isoform = transcript_id == focal
        decision = decision_by_id.get(internal_id, "review")
        fill, stroke = STATUS[decision]
        has_signal = signal.get(internal_id)
        if has_signal is False:
            stroke = NO_SIGNAL
        if not focal_isoform:
            fill, stroke = pale(fill), pale(stroke, 0.48)

        bases = sum(record["end"] - record["start"] + 1 for record in segments)
        state = "signal" if has_signal is True else "no signal" if has_signal is False else "signal unknown"
        svg += [
            f'<text x="{left - 8}" y="{y + 15}" text-anchor="end" class="meta">{html.escape(internal_id + (" CURRENT" if focal_isoform else ""))}</text>',
            f'<text x="{left + draw + 8}" y="{y + 15}" class="meta">{bases:,} bp · {bases // 3:,} aa · {html.escape(decision)} · {state}</text>',
        ]
        centre_y = y + cds_height / 2
        for first, second in zip(segments, segments[1:]):
            svg.append(
                f'<line x1="{x(first["end"] + 1):.2f}" y1="{centre_y:.2f}" x2="{x(second["start"]):.2f}" y2="{centre_y:.2f}" stroke="{stroke}" stroke-width="{2.4 if focal_isoform else 1}"/>'
            )
        for segment in segments:
            svg.append(
                block(
                    x(segment["start"]),
                    x(segment["end"] + 1),
                    y,
                    cds_height,
                    transcript["strand"],
                    fill,
                    stroke,
                    2.4 if focal_isoform else 1,
                )
            )

    svg += [
        f'<text x="{left}" y="{height - 8}" class="sub">Orange CDS outlines indicate has_signal=FALSE.</text>',
        "</svg>",
    ]
    out_path.write_text("\n".join(svg) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff3", type=Path)
    parser.add_argument("output_svg", type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--decisions", required=True, type=Path)
    parser.add_argument("--focal-internal-id", required=True)
    args = parser.parse_args()
    render(args.gff3, args.features, args.decisions, args.output_svg, args.focal_internal_id)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
