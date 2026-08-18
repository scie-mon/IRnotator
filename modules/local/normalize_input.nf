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


process NORMALIZE_MULTI_GENOME_GFF {
    tag "normalize-multi-genome-gff"

    container params.normalizer_container

    publishDir "${params.outdir}/input", mode: 'copy'

    input:
    path genome_fastas
    path annotation_gffs
    path normalizer_script
    path multi_normalizer_script

    output:
    tuple path("normalized_proteins.faa"),
          path("sequence_registry.tsv"),
          emit: normalized
    path "translation_report.tsv", emit: translation_report
    path "duplicate_provenance.tsv", emit: duplicate_provenance

    script:
    """
    python ${multi_normalizer_script} \
        --genome-fasta ${genome_fastas.join(' ')} \
        --annot-gff ${annotation_gffs.join(' ')} \
        --normalizer ${normalizer_script} \
        --out-faa normalized_proteins.faa \
        --out-registry sequence_registry.tsv \
        --out-report translation_report.tsv \
        --out-duplicate-provenance duplicate_provenance.tsv \
        --gff-protein-attribute '${params.gff_protein_attribute}' \
        --translation-table ${params.translation_table} \
        ${params.allow_internal_stops ? '--allow-internal-stops' : ''}
    """
}
