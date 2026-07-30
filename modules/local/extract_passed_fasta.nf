process EXTRACT_PASSED_FASTA {
    tag "extract_passed_fasta"

    publishDir "${params.outdir}/classification", mode: 'copy'

    input:
    path decision_tsv
    path hmm_hit_faa

    output:
    path "passed.faa", emit: fasta

    script:
    """
    python ${projectDir}/bin/extract_passed_fasta.py \\
        ${decision_tsv} \\
        ${hmm_hit_faa} \\
        passed.faa
    """
}
