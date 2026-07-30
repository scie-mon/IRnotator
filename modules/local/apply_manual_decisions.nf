process APPLY_MANUAL_DECISIONS {
    tag "apply_manual_decisions"

    publishDir "${params.outdir}/classification", mode: 'copy'

    input:
    path seed_decision_tsv
    path review_tsv

    output:
    path "decision.tsv", emit: decisions

    script:
    """
    python ${projectDir}/bin/apply_manual_decisions.py \\
        ${seed_decision_tsv} \\
        ${review_tsv} \\
        decision.tsv
    """
}
