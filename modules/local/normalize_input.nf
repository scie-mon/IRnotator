process NORMALIZE_PROTEINS {
    tag "normalize: ${proteins_faa.baseName}"
    container params.normalizer_container

    publishDir "${params.outdir}/input", mode: 'copy'

    input:
    path proteins_faa
    path normalizer_script

    output:
    tuple path("normalized_proteins.faa"), path("sequence_registry.tsv"), emit: normalized

    script:
    """
    python normalize_input.py \
        --proteins-faa ${proteins_faa} \
        --out-faa normalized_proteins.faa \
        --registry sequence_registry.tsv
    """
}


process NORMALIZE_GENOME_GFF {
    tag "normalize: ${genome_fasta.baseName}"
    container params.normalizer_container

    publishDir "${params.outdir}/input", mode: 'copy'

    input:
    path genome_fasta
    path annot_gff
    path normalizer_script

    output:
    tuple path("normalized_proteins.faa"), path("sequence_registry.tsv"), emit: normalized

    script:
    """
    python normalize_input.py \
        --genome-fasta ${genome_fasta} \
        --annot-gff ${annot_gff} \
        --out-faa normalized_proteins.faa \
        --registry sequence_registry.tsv \
        --gff-protein-attribute '${params.gff_protein_attribute}' \
        --translation-table ${params.translation_table} \
        ${params.allow_internal_stops ? '--allow-internal-stops' : ''}
    """
}


process NORMALIZE_PROTEINS_GFF {
    tag "normalize: ${proteins_faa.baseName}"
    container params.normalizer_container

    publishDir "${params.outdir}/input", mode: 'copy'

    input:
    path proteins_faa
    path annot_gff
    path normalizer_script

    output:
    tuple path("normalized_proteins.faa"), path("sequence_registry.tsv"), emit: normalized

    script:
    """
    python normalize_input.py \
        --proteins-faa ${proteins_faa} \
        --annot-gff ${annot_gff} \
        --out-faa normalized_proteins.faa \
        --registry sequence_registry.tsv \
        --gff-protein-attribute '${params.gff_protein_attribute}'
    """
}
