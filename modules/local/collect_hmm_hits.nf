process COLLECT_HMM_HITS {
    tag "collect_hmm_hits"

    publishDir "${params.outdir}/hmm_hits", mode: 'copy'

    input:
    path tbl_files

    output:
    path "HMM_ids.txt"

    script:
    """
    awk '!/^#/ && NF > 0 {print \$1}' ${tbl_files} | sort -u > HMM_ids.txt
    """
}
