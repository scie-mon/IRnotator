process EXTRACT_FASTA_BY_ID {
    tag "${fasta_in.baseName}"

    publishDir "${params.outdir}/hmm_hits", mode: 'copy'

    input:
    path ids_file
    path fasta_in

    output:
    path "HMM_hits.aa"

    script:
    """
    python ${projectDir}/bin/extract_fasta_by_id.py ${ids_file} ${fasta_in} HMM_hits.aa
    """
}
