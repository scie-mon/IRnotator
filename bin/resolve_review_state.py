#!/usr/bin/env python3
"""Resolve immutable merged-GFF revision state before a Nextflow run."""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

REVISION_RE = re.compile(r"^merged\.rev(\d+)\.gff3$")


def fail(message: str) -> None:
    raise SystemExit(f"resolve_review_state: {message}")


def revision_files(merged_dir: Path) -> dict[int, Path]:
    found: dict[int, Path] = {}
    for path in merged_dir.iterdir():
        match = REVISION_RE.match(path.name)
        if not match:
            continue
        if not path.is_file():
            fail(f"revision path is not a regular file: {path}")
        revision = int(match.group(1))
        if revision in found:
            fail(f"duplicate revision {revision}: {path}")
        found[revision] = path.resolve()
    return found


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--review-tsv", type=Path)
    args = parser.parse_args()

    merged_dir = args.merged_dir.expanduser().resolve()
    merged_dir.mkdir(parents=True, exist_ok=True)
    if not merged_dir.is_dir():
        fail(f"not a directory: {merged_dir}")

    unfiltered_gff = merged_dir / "merged.unfiltered.gff3"
    revisions = revision_files(merged_dir)
    revision_numbers = sorted(revisions)

    if revision_numbers and not unfiltered_gff.is_file():
        fail(
            "found merged.rev*.gff3 files but missing required "
            f"unfiltered reference: {unfiltered_gff}"
        )

    if revision_numbers:
        expected = list(range(revision_numbers[-1] + 1))
        if revision_numbers != expected:
            fail(
                "revision numbers must be contiguous from rev0; found "
                + ", ".join(f"rev{n}" for n in revision_numbers)
            )

    if args.review_tsv is not None:
        review_tsv = args.review_tsv.expanduser().resolve()
        if not review_tsv.is_file():
            fail(f"review TSV is not a regular file: {review_tsv}")
    else:
        review_tsv = None

    if not revision_numbers:
        mode = "bootstrap"
        revision_to_write = 0
        reviewer_context_kind = "current_unfiltered"
        reviewer_gff = unfiltered_gff.resolve()
    elif review_tsv is not None:
        mode = "advance"
        revision_to_write = revision_numbers[-1] + 1
        reviewer_context_kind = "previous_revision"
        reviewer_gff = revisions[revision_numbers[-1]]
    else:
        mode = "reuse"
        revision_to_write = None
        reviewer_context_kind = "existing_revision"
        reviewer_gff = revisions[revision_numbers[-1]]

    target_filtered_gff = (
        merged_dir / f"merged.rev{revision_to_write}.gff3"
        if revision_to_write is not None
        else None
    )
    if target_filtered_gff is not None and target_filtered_gff.exists():
        fail(f"refusing to overwrite existing revision: {target_filtered_gff}")

    state = {
        "mode": mode,
        "merged_dir": str(merged_dir),
        "unfiltered_gff": str(unfiltered_gff.resolve()),
        "review_tsv": str(review_tsv) if review_tsv else None,
        "existing_max_revision": revision_numbers[-1] if revision_numbers else None,
        "reviewer_context_kind": reviewer_context_kind,
        "reviewer_gff": str(reviewer_gff),
        "revision_to_write": revision_to_write,
        "target_filtered_gff": str(target_filtered_gff.resolve()) if target_filtered_gff else None,
    }

    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    os.replace(temporary, output)
    print(json.dumps(state, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as exc:
        print(f"resolve_review_state: {exc}", file=sys.stderr)
        raise SystemExit(1)
