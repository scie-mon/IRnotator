process SEED_DECISIONS {
    tag "seed_decisions"

    publishDir "${params.outdir}/classification", mode: 'copy',
        pattern: 'candidate_summary.tsv'
    // decision.tsv: channel only; final file published by apply (or copy step when no review)

    input:
    path features_tsv
    path hardcoded_tsv
    val  manual_review

    output:
    path "candidate_summary.tsv", emit: summary
    path "decision.tsv",          emit: decisions

    script:
    """
    python ${projectDir}/bin/seed_decisions.py \\
        ${features_tsv} \\
        ${hardcoded_tsv} \\
        candidate_summary.tsv \\
        decision.tsv \\
        ${manual_review}
    """
}
