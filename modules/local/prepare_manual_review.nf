process PREPARE_MANUAL_REVIEW {
    tag "prepare_manual_review:${reviewer_revision}:${written_revision}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path deeptmhmm_dirs
    path scores_tsv
    path registry_tsv
    path features_tsv
    path effective_decisions
    path applied_review_tsv
    path reviewer_gff
    path reviewer_html
    path unavailable_context_svg
    path renderer_script
    val reviewer_revision
    val written_revision
    val genome_context_available

    output:
    path "manual_review", emit: review_dir
    path "manual_review/review-manifest.json", emit: manifest
    path "manual_review/topology-reviewer.html", emit: reviewer

    script:
    """
    set -euo pipefail

    mkdir -p manual_review/plots manual_review/genes
    cp ${reviewer_html} manual_review/topology-reviewer.html
    cp ${scores_tsv} manual_review/hardcoded_scores.tsv
    cp ${registry_tsv} manual_review/sequence_registry.tsv
    cp ${features_tsv} manual_review/deeptmhmm_features.tsv
    cp ${effective_decisions} manual_review/effective_decisions.tsv

    if [[ ${genome_context_available} == true ]]; then
        cp ${reviewer_gff} manual_review/reviewer_context.gff3
    else
        printf '%s\\n' '## genome-context-unavailable: protein-only input' > manual_review/reviewer_context.gff3
        cp ${unavailable_context_svg} manual_review/genes/genome-context-unavailable.svg
    fi

    python3 - <<'PY'
    from __future__ import annotations

    import csv
    import json
    import shutil
    import subprocess
    from pathlib import Path

    REVIEW_DIR = Path("manual_review")
    VALID_DECISIONS = {"accept", "reject", "review"}
    GENOME_CONTEXT_AVAILABLE = "${genome_context_available}".lower() == "true"


    def fail(message: str) -> None:
        raise SystemExit(f"PREPARE_MANUAL_REVIEW: {message}")


    def read_tsv(path: Path) -> list[dict[str, str]]:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            if not reader.fieldnames:
                fail(f"empty TSV: {path}")
            return list(reader)


    def attrs(text: str) -> dict[str, str]:
        return {
            key: value
            for item in text.rstrip(";").split(";")
            if "=" in item
            for key, value in [item.split("=", 1)]
        }


    def parse_context(path: Path) -> dict[str, dict[str, object]]:
        items = {}
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\\n").split("\\t")
                if len(fields) != 9:
                    fail(f"{path}:{line_number}: expected 9-column GFF3")
                contig, _, kind, start, end, _, strand, _, raw_attrs = fields
                if kind not in {"mRNA", "transcript"}:
                    continue
                parsed = attrs(raw_attrs)
                internal_id = parsed.get("internal_id", "").strip()
                transcript_id = parsed.get("ID", "").strip()
                parents = [parent for parent in parsed.get("Parent", "").split(",") if parent]
                if not internal_id or not transcript_id or len(parents) != 1:
                    fail(f"{path}:{line_number}: mRNA requires internal_id, ID, and one Parent")
                if internal_id in items:
                    fail(f"duplicate internal_id in reviewer context GFF: {internal_id}")
                items[internal_id] = {
                    "seq_id": internal_id,
                    "gene_id": parents[0],
                    "transcript_id": transcript_id,
                    "contig": contig,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                }
        if not items:
            fail(f"no mRNA records found in reviewer context GFF: {path}")
        return items


    def normalize_decision(decision: str) -> str:
        return "undecided" if decision == "review" else decision


    registry_by_id = {}
    for row in read_tsv(REVIEW_DIR / "sequence_registry.tsv"):
        internal_id = row.get("internal_id", "").strip()
        if not internal_id or internal_id in registry_by_id:
            fail(f"invalid or duplicate internal_id in sequence registry: {internal_id!r}")
        registry_by_id[internal_id] = row

    reason_by_id = {}
    for row in read_tsv(REVIEW_DIR / "hardcoded_scores.tsv"):
        internal_id = row.get("protein_id", "").strip()
        if internal_id:
            reason_by_id[internal_id] = row.get("hardcoded_reason", "").strip()

    decision_by_id = {}
    for row in read_tsv(REVIEW_DIR / "effective_decisions.tsv"):
        internal_id = row.get("protein_id", "").strip()
        decision = row.get("final_decision", "").strip().lower()
        if not internal_id:
            fail("empty protein_id in effective decisions TSV")
        if decision not in VALID_DECISIONS:
            fail(f"invalid final_decision for {internal_id}: {decision!r}")
        if internal_id in decision_by_id:
            fail(f"duplicate protein_id in effective decisions TSV: {internal_id}")
        decision_by_id[internal_id] = decision

    if GENOME_CONTEXT_AVAILABLE:
        context_by_id = parse_context(REVIEW_DIR / "reviewer_context.gff3")
    else:
        context_by_id = {}
        for internal_id in decision_by_id:
            registry = registry_by_id.get(internal_id)
            if registry is None:
                fail(f"effective-decision internal_id is absent from registry: {internal_id}")
            source_id = registry.get("source_id", "").strip()
            source_header = registry.get("source_header", "").strip()
            context_by_id[internal_id] = {
                "seq_id": internal_id,
                "gene_id": source_id or internal_id,
                "transcript_id": source_header or internal_id,
                "contig": "protein-only",
                "start": 0,
                "end": 0,
                "strand": ".",
            }

    items = []
    for internal_id, item in context_by_id.items():
        if internal_id not in registry_by_id:
            fail(f"review-context internal_id is absent from registry: {internal_id}")
        if internal_id not in decision_by_id:
            fail(f"review-context internal_id is absent from effective decisions: {internal_id}")

        plots = list(Path(internal_id).rglob("*.png"))
        if len(plots) != 1:
            fail(f"expected one DeepTMHMM PNG for {internal_id}; found {len(plots)}")

        plot_out = REVIEW_DIR / "plots" / f"{internal_id}.png"
        shutil.copy2(plots[0], plot_out)

        if GENOME_CONTEXT_AVAILABLE:
            gene_out = REVIEW_DIR / "genes" / f"{internal_id}.svg"
            subprocess.run([
                "python3", "${renderer_script}", "${reviewer_gff}", str(gene_out),
                "--features", "${features_tsv}",
                "--decisions", "${effective_decisions}",
                "--focal-internal-id", internal_id,
            ], check=True)
            gene_image = f"genes/{internal_id}.svg"
        else:
            gene_image = "genes/genome-context-unavailable.svg"

        registry = registry_by_id[internal_id]
        item["seed"] = normalize_decision(decision_by_id[internal_id])
        item.update({
            "source_id": registry.get("source_id", ""),
            "source_header": registry.get("source_header", ""),
            "reason": reason_by_id.get(internal_id, ""),
            "topology_image": f"plots/{internal_id}.png",
            "gene_image": gene_image,
        })
        items.append(item)

    items.sort(key=lambda item: (
        item["contig"], item["start"], item["gene_id"], item["end"], item["transcript_id"],
    ))
    for order, item in enumerate(items, 1):
        item["review_order"] = order

    (REVIEW_DIR / "review-manifest.json").write_text(json.dumps({
        "reviewer_revision": ${reviewer_revision},
        "genome_context_available": GENOME_CONTEXT_AVAILABLE,
        "items": items,
    }, indent=2) + "\\n")
    PY

    if [[ ${written_revision} -gt 0 ]]; then
        cp ${applied_review_tsv} manual_review/review.applied.rev${written_revision}.tsv
        rm -f "${params.outdir}/manual_review/review.tsv"
    fi

    cat > manual_review/rerun_after_review.sh <<'EOF'
    #!/usr/bin/env bash
    set -euo pipefail

    # Remove cached tasks that consume review state or generate review assets.
    # Usage: ./rerun_after_review.sh [WORK_DIR]
    work_dir="\${1:-work}"
    if [[ ! -d "\${work_dir}" ]]; then
        echo "Work directory not found: \${work_dir}" >&2
        exit 1
    fi

    remove_cached_tasks() {
        local task_name="\${1}"
        find "\${work_dir}" -type f -name "\${task_name}" -print0 |
        while IFS= read -r -d '' task_file; do
            if grep -qE 'merge_ir_isoforms.py|render_gene_cds.py|reviewer_context.gff3|review.tsv' "\${task_file}"; then
                task_dir="\$(dirname "\${task_file}")"
                echo "Removing cached review-dependent task: \${task_dir}"
                rm -rf "\${task_dir}"
            fi
        done
    }

    remove_cached_tasks .command.sh
    remove_cached_tasks .command.run
    EOF
    chmod +x manual_review/rerun_after_review.sh

    cat > manual_review/README.txt <<'EOF'
    Manual review package
    =====================

    The manifest is ordered by merged gene model, then isoform.

    1. Serve this directory: python3 -m http.server 8000
    2. Open topology-reviewer.html.
    3. Export review.tsv into this directory.
    4. Run ./rerun_after_review.sh /path/to/IRnotator/work
    5. Re-run the same Nextflow command with -resume.

    seq_id is the IRnotator internal ID and is the only decision key.
    EOF
    """.stripIndent()
}
