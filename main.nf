nextflow.enable.dsl = 2

include { NORMALIZE_PROTEINS }         from './modules/local/normalize_input'
include { NORMALIZE_GENOME_GFF }       from './modules/local/normalize_input'
include { NORMALIZE_PROTEINS_GFF }     from './modules/local/normalize_input'

include { HMMSEARCH_IR }               from './modules/local/hmmsearch_ir'
include { COLLECT_HMM_HITS }           from './modules/local/collect_hmm_hits'
include { EXTRACT_FASTA_BY_ID }        from './modules/local/extract_fasta_by_id'
include { SPLIT_FASTA_BY_SEQID }       from './modules/local/split_fasta_by_seqid'
include { DEEPTMHMM_TOPOLOGY }         from './modules/local/deeptmhmm_topology'
include { COLLECT_DEEPTMHMM_SUMMARY }  from './modules/local/collect_deeptmhmm_summary'
include { PARSE_DEEPTMHMM_FEATURES }   from './modules/local/parse_deeptmhmm_features'
include { HARDCODED_CLASSIFIER }       from './modules/local/hardcoded_classifier'
include { SEED_DECISIONS }             from './modules/local/seed_decisions'
include { PREPARE_MANUAL_REVIEW }      from './modules/local/prepare_manual_review'
include { GATE_REVIEW_TSV }            from './modules/local/gate_review_tsv'
include { APPLY_MANUAL_DECISIONS }     from './modules/local/apply_manual_decisions'
include { EXTRACT_PASSED_FASTA }       from './modules/local/extract_passed_fasta'

params.proteins_faa          = params.proteins_faa ?: null
params.genome_fasta          = params.genome_fasta ?: null
params.annot_gff             = params.annot_gff ?: null
params.gff_protein_attribute = params.gff_protein_attribute ?: 'protein_id'
params.translation_table     = params.translation_table ?: 1
params.allow_internal_stops  = params.allow_internal_stops ?: false

params.hmm_dir       = params.hmm_dir ?: 'test/hmms'
params.outdir        = params.outdir ?: 'results'
params.deeptmhmm_dir = params.deeptmhmm_dir ?: null
params.manual_review = params.manual_review ?: true

workflow {
    def has_proteins = params.proteins_faa as boolean
    def has_genome   = params.genome_fasta as boolean
    def has_gff      = params.annot_gff as boolean

    if (has_proteins && has_genome && has_gff) {
        throw new IllegalArgumentException(
            'Invalid input: do not provide proteins_faa, genome_fasta, and annot_gff together.'
        )
    }

    normalizer_script_ch = Channel.value(
        file("${projectDir}/bin/normalize_input.py")
    )

    if (has_proteins && !has_genome && !has_gff) {
        proteins_input_ch = Channel.fromPath(
            params.proteins_faa,
            checkIfExists: true
        )

        normalize_res = NORMALIZE_PROTEINS(
            proteins_input_ch,
            normalizer_script_ch
        )
    }
    else if (!has_proteins && has_genome && has_gff) {
        genome_input_ch = Channel.fromPath(
            params.genome_fasta,
            checkIfExists: true
        )
        gff_input_ch = Channel.fromPath(
            params.annot_gff,
            checkIfExists: true
        )

        normalize_res = NORMALIZE_GENOME_GFF(
            genome_input_ch,
            gff_input_ch,
            normalizer_script_ch
        )
    }
    else if (has_proteins && !has_genome && has_gff) {
        proteins_input_ch = Channel.fromPath(
            params.proteins_faa,
            checkIfExists: true
        )
        gff_input_ch = Channel.fromPath(
            params.annot_gff,
            checkIfExists: true
        )

        normalize_res = NORMALIZE_PROTEINS_GFF(
            proteins_input_ch,
            gff_input_ch,
            normalizer_script_ch
        )
    }
    else {
        throw new IllegalArgumentException(
            'Invalid input: provide exactly one valid combination: ' +
            '(1) --proteins_faa; ' +
            '(2) --genome_fasta plus --annot_gff; or ' +
            '(3) --proteins_faa plus --annot_gff.'
        )
    }

    /*
     * Canonical identity boundary:
     * normalized_proteins.faa has only IRN_* identifiers.
     * All downstream analysis uses IRN_* exclusively.
     * sequence_registry.tsv retains source provenance and output naming metadata.
     */
    normalized_proteins_ch = normalize_res.normalized.map { proteins_faa, registry_tsv ->
        proteins_faa
    }

    sequence_registry_ch = normalize_res.normalized.map { proteins_faa, registry_tsv ->
        registry_tsv
    }

    hmms_ch          = Channel.fromPath("${params.hmm_dir}/*.hmm", checkIfExists: true)
    deeptmhmm_dir_ch = Channel.value(file(params.deeptmhmm_dir))

    hmm_input_ch = hmms_ch.combine(normalized_proteins_ch)

    hmm_res     = HMMSEARCH_IR(hmm_input_ch)
    hmm_tbls_ch = hmm_res.tbl

    hmm_ids_ch     = COLLECT_HMM_HITS(hmm_tbls_ch.collect())
    hmm_hit_faa_ch = EXTRACT_FASTA_BY_ID(hmm_ids_ch, normalized_proteins_ch)

    split_res         = SPLIT_FASTA_BY_SEQID(hmm_hit_faa_ch)
    single_seq_faa_ch = split_res.fasta_files.flatten()

    deeptmhmm_res = DEEPTMHMM_TOPOLOGY(
        single_seq_faa_ch,
        deeptmhmm_dir_ch
    )

    summary_res = COLLECT_DEEPTMHMM_SUMMARY(
        deeptmhmm_res.results.collect()
    )

    parsed_res = PARSE_DEEPTMHMM_FEATURES(
        summary_res.summary_dir
    )

    hardcoded_res = HARDCODED_CLASSIFIER(parsed_res.features)

    seed_res = SEED_DECISIONS(
        parsed_res.features,
        hardcoded_res.scores,
        params.manual_review
    )

    if (params.manual_review) {
        prepare_res = PREPARE_MANUAL_REVIEW(
            deeptmhmm_res.results.collect(),
            hardcoded_res.scores,
            sequence_registry_ch,
            file("${projectDir}/assets/reviewer/topology-reviewer.html")
        )

        review_tsv_path = file(
            "${params.outdir}/manual_review/review.tsv"
        ).toAbsolutePath().toString()

        gate_res = GATE_REVIEW_TSV(
            review_tsv_path,
            prepare_res.review_dir
        )

        decisions_ch = APPLY_MANUAL_DECISIONS(
            seed_res.decisions,
            gate_res.review_tsv
        ).decisions
    }
    else {
        decisions_ch = seed_res.decisions
    }

    passed_res = EXTRACT_PASSED_FASTA(
        decisions_ch,
        hmm_hit_faa_ch
    )

    // Future: join passed_res.fasta to sequence_registry_ch,
    // assign final_output_id, and create final FASTA/GFF/report outputs.
}
