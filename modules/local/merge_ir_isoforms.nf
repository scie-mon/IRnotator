process MERGE_IR_ISOFORMS {
    tag "merge_ir_isoforms"

    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    path annotation_gffs
    path sequence_registry
    path decisions
    path merger_script

    output:
    path "merged_IR.gff3", emit: merged_gff
    path "merged_isoform_audit.tsv", emit: audit

    script:
    """
    python ${merger_script} \
        --gff ${annotation_gffs.join(' ')} \
        --registry ${sequence_registry} \
        --decisions ${decisions} \
        --overlap-fraction ${params.isoform_overlap_fraction} \
        --out-gff merged_IR.gff3 \
        --out-audit merged_isoform_audit.tsv
    """
}
