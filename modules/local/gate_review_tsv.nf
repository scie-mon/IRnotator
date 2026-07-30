process GATE_REVIEW_TSV {
    tag "gate_review_tsv"

    input:
    val review_tsv_path   // absolute/resolved path string under outdir
    val _prepare_done     // dependency: PREPARE_MANUAL_REVIEW finished

    output:
    path "review.tsv", emit: review_tsv

    script:
    """
    if [ ! -f "${review_tsv_path}" ]; then
        echo "ERROR: Manual review file not found:" >&2
        echo "  ${review_tsv_path}" >&2
        echo "" >&2
        echo "Export review.tsv from the topology reviewer (key e)," >&2
        echo "place it at that path, then re-run with -resume." >&2
        exit 1
    fi
    cp "${review_tsv_path}" review.tsv
    """
}
