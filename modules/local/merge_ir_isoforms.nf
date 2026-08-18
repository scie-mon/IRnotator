process MERGE_IR_ISOFORMS {
    tag "${merge_mode}:${gff_name}"

    publishDir "${params.outdir}/merged_gff", mode: 'copy', overwrite: false

    input:
    path annotation_gffs
    path sequence_registry
    path decisions
    path duplicate_provenance
    path merger_script
    val merge_mode
    val gff_name
    val audit_name

    output:
    path "${gff_name}", emit: merged_gff
    path "${audit_name}", emit: audit

    script:
    """
    set -euo pipefail

    python3 ${merger_script} \\
        --gff ${annotation_gffs.join(' ')} \\
        --registry ${sequence_registry} \\
        --decisions ${decisions} \\
        --duplicate-provenance ${duplicate_provenance} \\
        --overlap-fraction ${params.isoform_overlap_fraction} \\
        --mode ${merge_mode} \\
        --out-gff ${gff_name} \\
        --out-audit ${audit_name}
    """
}
