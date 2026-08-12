#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from urllib.parse import unquote


def attrs(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[unquote(key)] = unquote(value)
    return result


def track_name(path: Path) -> str:
    return re.sub(r"\.(?:gff3?|gtf)(?:\.gz)?$", "", path.name, flags=re.I)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def span(model: dict) -> tuple[int, int]:
    return min(item["start"] for item in model["cds"]), max(item["end"] for item in model["cds"])


def overlap(a: tuple[int, int], b: tuple[int, int]) -> float:
    left, right = max(a[0], b[0]), min(a[1], b[1])
    if right < left:
        return 0.0
    return (right - left + 1) / min(a[1] - a[0] + 1, b[1] - b[0] + 1)


def parse_gffs(paths: list[Path]) -> dict[tuple[str, str], dict]:
    models = {}
    for path in paths:
        track = track_name(path)
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9 or fields[2] != "CDS":
                    continue
                contig, source, _, start, end, score, strand, phase, raw = fields
                parsed = attrs(raw)
                parents = [item for item in parsed.get("Parent", "").split(",") if item]
                if len(parents) != 1:
                    raise ValueError(f"{path}:{line_number}: CDS must have one Parent")
                key = (track, parents[0])
                model = models.setdefault(key, {
                    "track": track,
                    "source": source,
                    "transcript_id": parents[0],
                    "contig": contig,
                    "strand": strand,
                    "cds": [],
                })
                if (model["contig"], model["strand"]) != (contig, strand):
                    raise ValueError(f"inconsistent CDS location for {key}")
                model["cds"].append({
                    "start": int(start),
                    "end": int(end),
                    "score": score,
                    "strand": strand,
                    "phase": phase,
                    "attrs": parsed,
                })
    return models


def clusters(candidates: list[dict], threshold: float) -> list[list[int]]:
    result = []
    remaining = set(range(len(candidates)))
    while remaining:
        group = {remaining.pop()}
        changed = True
        while changed:
            changed = False
            for index in list(remaining):
                if any(
                    candidates[index]["contig"] == candidates[other]["contig"]
                    and candidates[index]["strand"] == candidates[other]["strand"]
                    and overlap(span(candidates[index]), span(candidates[other])) >= threshold
                    for other in group
                ):
                    group.add(index)
                    remaining.remove(index)
                    changed = True
        result.append(sorted(group))
    return result


def source_attributes(models: list[dict]) -> str:
    values = set()
    for model in models:
        for segment in model["cds"]:
            for key, value in segment["attrs"].items():
                if key != "translation":
                    values.add(
                        f"{key}={value}"
                        .replace(";", "%3B")
                        .replace("=", "%3D")
                        .replace(",", "%2C")
                    )
    return "%2C".join(sorted(values))


def candidates_from_registry(
    registry: list[dict[str, str]],
    decision_by_id: dict[str, dict[str, str]],
    models: dict[tuple[str, str], dict],
    mode: str,
) -> list[dict]:
    candidates = []
    for row in registry:
        internal_id = row["internal_id"].strip()
        decision = decision_by_id.get(internal_id)
        if decision is None:
            continue
        final_decision = decision["final_decision"].strip().lower()
        if final_decision not in {"accept", "reject", "review"}:
            raise ValueError(f"invalid final_decision for {internal_id}: {final_decision!r}")
        if mode == "filtered" and final_decision != "accept":
            continue
        key = (row["annotation_track"], row["source_id"])
        if key not in models:
            raise ValueError(f"no GFF CDS model for {internal_id}: {key}")
        candidate = dict(models[key])
        candidate.update(
            internal_id=internal_id,
            topology_class=decision["topology_class"].strip(),
            final_decision=final_decision,
            registry=row,
        )
        candidates.append(candidate)
    return candidates


def select_models(candidates: list[dict], threshold: float, mode: str) -> list[list[dict]]:
    initial_groups = clusters(candidates, threshold)
    if mode == "unfiltered":
        return [[candidates[index] for index in group] for group in initial_groups]

    selected = []
    for group in initial_groups:
        strict = [candidates[index] for index in group if candidates[index]["topology_class"] == "strict"]
        selected.extend(strict if strict else [candidates[index] for index in group])

    if not selected:
        return []
    return [[selected[index] for index in group] for group in clusters(selected, threshold)]


def write_outputs(groups: list[list[dict]], out_gff: Path, out_audit: Path) -> None:
    audit = []
    with out_gff.open("w") as out:
        out.write("##gff-version 3\n")
        for models_out in groups:
            first = models_out[0]
            contig, strand = first["contig"], first["strand"]
            gene_start = min(span(model)[0] for model in models_out)
            gene_end = max(span(model)[1] for model in models_out)
            gene_id = f"{contig}_{gene_start}-{gene_end}{'f' if strand == '+' else 'r'}"
            source_attr = source_attributes(models_out)
            out.write(
                f"{contig}\tMerged\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
                f"ID={gene_id};Name={gene_id};source_attributes={source_attr}\n"
            )
            ordered = sorted(models_out, key=lambda model: (*span(model), model["internal_id"]))
            multiple = len(ordered) > 1
            for iso_index, model in enumerate(ordered, 1):
                start, end = span(model)
                label = f"{gene_id}_X{iso_index}" if multiple else gene_id
                mrna_id = f"{label}-{start}-{end}-{model['internal_id']}"
                source = model["source"] or model["track"]
                out.write(
                    f"{contig}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
                    f"ID={mrna_id};Parent={gene_id};tm={model['topology_class']};"
                    f"decision={model['final_decision']};source={source};"
                    f"source_attributes={source_attr};Name={label} mRNA;"
                    f"internal_id={model['internal_id']};annotation_track={model['track']}\n"
                )
                for segment in sorted(model["cds"], key=lambda item: item["start"]):
                    translation = segment["attrs"].get("translation", "")
                    extra = f";translation={translation}" if translation else ""
                    out.write(
                        f"{contig}\t{source}\tCDS\t{segment['start']}\t{segment['end']}\t.\t"
                        f"{strand}\t{segment['phase']}\tID={mrna_id}.CDS;Parent={mrna_id};"
                        f"Name={label} CDS{extra};source_attributes={source_attr};"
                        f"internal_id={model['internal_id']}\n"
                    )
                audit.append({
                    "merged_gene_id": gene_id,
                    "merged_isoform_id": mrna_id,
                    "internal_id": model["internal_id"],
                    "final_decision": model["final_decision"],
                    "annotation_track": model["track"],
                    "source_id": model["registry"]["source_id"],
                    "transcript_id": model["registry"].get("transcript_id", ""),
                    "topology_class": model["topology_class"],
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                })

    fields = [
        "merged_gene_id", "merged_isoform_id", "internal_id", "final_decision",
        "annotation_track", "source_id", "transcript_id", "topology_class",
        "contig", "start", "end", "strand",
    ]
    with out_audit.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", nargs="+", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--decisions", required=True)
    parser.add_argument("--overlap-fraction", type=float, default=0.60)
    parser.add_argument("--mode", choices=("unfiltered", "filtered"), default="filtered")
    parser.add_argument("--out-gff", required=True)
    parser.add_argument("--out-audit", required=True)
    args = parser.parse_args()

    if not 0 <= args.overlap_fraction <= 1:
        raise ValueError("--overlap-fraction must be between 0 and 1")

    registry = read_tsv(Path(args.registry))
    decisions = read_tsv(Path(args.decisions))
    decision_by_id = {row["protein_id"].strip(): row for row in decisions}
    models = parse_gffs([Path(path) for path in args.gff])
    candidates = candidates_from_registry(registry, decision_by_id, models, args.mode)
    groups = select_models(candidates, args.overlap_fraction, args.mode)
    write_outputs(groups, Path(args.out_gff), Path(args.out_audit))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
