process SEED_DECISIONS {
    tag "seed_decisions"

    publishDir "${params.outdir}/classification", mode: 'copy'

    input:
    path features_tsv
    path hardcoded_tsv

    output:
    path "candidate_summary.tsv", emit: summary
    path "decision.tsv", emit: decisions

    script:
    """
    python ${projectDir}/bin/seed_decisions.py \
        ${features_tsv} \
        ${hardcoded_tsv} \
        candidate_summary.tsv \
        decision.tsv
    """
}