process HMMSEARCH_IR {
    tag "${hmm_file.baseName}"

    publishDir "${params.outdir}/hmmsearch", mode: 'copy'

    container params.hmmer_container

    input:
    tuple path(hmm_file), path(proteins_faa)

    output:
    path "${hmm_file.baseName}.tbl", emit: tbl
    path "${hmm_file.baseName}.txt", emit: txt

    script:
    """
    hmmsearch \
        -E ${params.hmm_evalue} \
        --tblout ${hmm_file.baseName}.tbl \
        ${hmm_file} \
        ${proteins_faa} \
        > ${hmm_file.baseName}.txt
    """
}
