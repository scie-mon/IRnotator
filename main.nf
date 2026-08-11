nextflow.enable.dsl = 2

include { NORMALIZE_PROTEINS }                from './modules/local/normalize_input'
include { NORMALIZE_MULTI_GENOME_GFF }        from './modules/local/normalize_input'

include { HMMSEARCH_IR }                      from './modules/local/hmmsearch_ir'
include { COLLECT_HMM_HITS }                  from './modules/local/collect_hmm_hits'
include { EXTRACT_FASTA_BY_ID }               from './modules/local/extract_fasta_by_id'
include { DEEPTMHMM_TOPOLOGY }                from './modules/local/deeptmhmm_topology'
include { INITIALIZE_DEEPTMHMM_RESULTS }      from './modules/local/initialize_deeptmhmm_results'
include { FINALIZE_DEEPTMHMM_RESULTS }        from './modules/local/finalize_deeptmhmm_results'
include { ASSERT_COMPLETE_DEEPTMHMM_RESULTS } from './modules/local/finalize_deeptmhmm_results'
include { COLLECT_DEEPTMHMM_SUMMARY }         from './modules/local/collect_deeptmhmm_summary'
include { PARSE_DEEPTMHMM_FEATURES }          from './modules/local/parse_deeptmhmm_features'
include { HARDCODED_CLASSIFIER }              from './modules/local/hardcoded_classifier'
include { SEED_DECISIONS }                    from './modules/local/seed_decisions'
include { PREPARE_MANUAL_REVIEW }             from './modules/local/prepare_manual_review'
include { GATE_REVIEW_TSV }                   from './modules/local/gate_review_tsv'
include { APPLY_MANUAL_DECISIONS }            from './modules/local/apply_manual_decisions'
include { EXTRACT_PASSED_FASTA }              from './modules/local/extract_passed_fasta'
include { MERGE_IR_ISOFORMS }                 from './modules/local/merge_ir_isoforms'

params.proteins_faa             = params.proteins_faa ?: null
params.genome_fasta             = params.genome_fasta ?: null
params.annot_gff                = params.annot_gff ?: null
params.gff_protein_attribute    = params.gff_protein_attribute ?: 'protein_id'
params.translation_table        = params.translation_table ?: 1
params.allow_internal_stops     = params.allow_internal_stops ?: false
params.isoform_overlap_fraction = params.isoform_overlap_fraction ?: 0.60

params.hmm_dir = params.hmm_dir ?: 'test/hmms'
params.outdir = params.outdir ?: 'results'
params.deeptmhmm_dir = params.deeptmhmm_dir ?: null
params.run_deeptmhmm = params.run_deeptmhmm == null ? true : params.run_deeptmhmm
params.deeptmhmm_salvage_paths = params.deeptmhmm_salvage_paths ?: []
params.manual_review = params.manual_review == null ? true : params.manual_review
params.review_tsv = params.review_tsv ?: null

workflow {
    def has_proteins = params.proteins_faa as boolean
    def has_genome   = params.genome_fasta as boolean
    def has_gff      = params.annot_gff as boolean

    if (has_proteins && has_gff) {
        throw new IllegalArgumentException(
            'The protein+GFF mode is deprecated and unsupported. ' +
            'Use --proteins_faa alone or --genome_fasta with --annot_gff.'
        )
    }

    if (has_proteins && has_genome) {
        throw new IllegalArgumentException(
            'Invalid input: do not provide proteins_faa together with genome_fasta.'
        )
    }

    if (!has_proteins && !has_genome) {
        throw new IllegalArgumentException(
            'Invalid input: provide --proteins_faa or --genome_fasta plus --annot_gff.'
        )
    }

    if (has_genome && !has_gff) {
        throw new IllegalArgumentException(
            'Genome FASTA input requires one or more --annot_gff files.'
        )
    }

    if (params.isoform_overlap_fraction < 0.0 || params.isoform_overlap_fraction > 1.0) {
        throw new IllegalArgumentException(
            '--isoform_overlap_fraction must be between 0 and 1.'
        )
    }

    normalizer_script_ch = Channel.value(
        file("${projectDir}/bin/normalize_input.py")
    )

    multi_normalizer_script_ch = Channel.value(
        file("${projectDir}/bin/multi_normalize.py")
    )

    merger_script_ch = Channel.value(
        file("${projectDir}/bin/merge_ir_isoforms.py")
    )

    annotation_gffs_ch = null

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
        genome_input_ch = Channel
            .fromPath(params.genome_fasta, checkIfExists: true)
            .collect()

        annotation_gffs_ch = Channel
            .fromPath(params.annot_gff, checkIfExists: true)
            .collect()

        normalize_res = NORMALIZE_MULTI_GENOME_GFF(
            genome_input_ch,
            annotation_gffs_ch,
            normalizer_script_ch,
            multi_normalizer_script_ch
        )
    }
    else {
        throw new IllegalArgumentException(
            'Invalid input: supported combinations are ' +
            '(1) --proteins_faa; or ' +
            '(2) --genome_fasta plus one or more --annot_gff files.'
        )
    }

    normalized_proteins_ch = normalize_res.normalized.map {
        proteins_faa, registry_tsv -> proteins_faa
    }

    sequence_registry_ch = normalize_res.normalized.map {
        proteins_faa, registry_tsv -> registry_tsv
    }

    hmms_ch = Channel.fromPath(
        "${params.hmm_dir}/*.hmm",
        checkIfExists: true
    )

    hmm_input_ch = hmms_ch.combine(normalized_proteins_ch)
    hmm_res = HMMSEARCH_IR(hmm_input_ch)
    hmm_tbls_ch = hmm_res.tbl

    hmm_ids_ch = COLLECT_HMM_HITS(hmm_tbls_ch.collect())
    hmm_hit_faa_ch = EXTRACT_FASTA_BY_ID(
        hmm_ids_ch,
        normalized_proteins_ch
    )

    def salvage_paths = params.deeptmhmm_salvage_paths instanceof Collection \
        ? params.deeptmhmm_salvage_paths.collect { it.toString() } \
        : params.deeptmhmm_salvage_paths
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }

    initial_deeptmhmm_res = INITIALIZE_DEEPTMHMM_RESULTS(
        hmm_hit_faa_ch,
        Channel.value(salvage_paths)
    )

    initial_manifest_ch = initial_deeptmhmm_res.manifest.map { manifest_tsv ->
        def rows = manifest_tsv
            .toFile()
            .readLines()
            .drop(1)
            .findAll { it }
            .collect { it.split('\\t', -1) }

        def salvaged = rows.count { fields -> fields[3] == 'salvaged' }

        log.info(
            "DeepTMHMM candidates: ${rows.size()}; " +
            "salvaged: ${salvaged}; " +
            "queued for new prediction: ${rows.size() - salvaged}"
        )

        manifest_tsv
    }

    def generated_results_ch

    if (params.run_deeptmhmm) {
        if (!params.deeptmhmm_dir) {
            throw new IllegalArgumentException(
                'Missing --deeptmhmm_dir while --run_deeptmhmm is enabled.'
            )
        }

        def unresolved_deeptmhmm_faa_ch = initial_deeptmhmm_res
            .unresolved_fasta
            .flatten()

        deeptmhmm_res = DEEPTMHMM_TOPOLOGY(
            unresolved_deeptmhmm_faa_ch,
            Channel.value(file(params.deeptmhmm_dir))
        )

        generated_results_ch = deeptmhmm_res.results
            .map { sequence_id, result_dir -> result_dir }
            .collect()
            .ifEmpty { [] }
    }
    else {
        generated_results_ch = Channel.value([])
    }

    final_deeptmhmm_res = FINALIZE_DEEPTMHMM_RESULTS(
        initial_manifest_ch,
        generated_results_ch
    )

    ASSERT_COMPLETE_DEEPTMHMM_RESULTS(
        final_deeptmhmm_res.unresolved_fasta
    )

    final_deeptmhmm_dirs_ch = final_deeptmhmm_res.result_dirs.flatten()

    summary_res = COLLECT_DEEPTMHMM_SUMMARY(
        final_deeptmhmm_res.manifest,
        final_deeptmhmm_res.results_dir
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
            final_deeptmhmm_dirs_ch.collect(),
            hardcoded_res.scores,
            sequence_registry_ch,
            file("${projectDir}/assets/reviewer/topology-reviewer.html")
        )

        def review_file = params.review_tsv \
            ? file(params.review_tsv as String) \
            : file("${params.outdir}/manual_review/review.tsv")

        if (params.review_tsv && !review_file.exists()) {
            throw new IllegalArgumentException(
                "Manual review TSV not found: ${review_file.toAbsolutePath()}"
            )
        }

        if (review_file.exists()) {
            gate_res = GATE_REVIEW_TSV(
                review_file.toAbsolutePath().toString(),
                prepare_res.review_dir
            )

            final_decisions_ch = APPLY_MANUAL_DECISIONS(
                seed_res.decisions,
                gate_res.review_tsv
            ).decisions

            EXTRACT_PASSED_FASTA(
                final_decisions_ch,
                hmm_hit_faa_ch
            )

            MERGE_IR_ISOFORMS(
                annotation_gffs_ch,
                sequence_registry_ch,
                final_decisions_ch,
                merger_script_ch
            )
        }
        else {
            log.info """
            Manual review package ready:
              ${params.outdir}/manual_review

            Next steps:
              1. cd ${params.outdir}/manual_review
              2. python3 -m http.server 8000
              3. Open http://localhost:8000/topology-reviewer.html
              4. Export review.tsv into this directory
              5. Re-run the same Nextflow command with -resume

            Optional:
              --review_tsv /absolute/path/to/review.tsv
            """.stripIndent()
        }
    }
    else {
        final_decisions_ch = seed_res.decisions

        EXTRACT_PASSED_FASTA(
            final_decisions_ch,
            hmm_hit_faa_ch
        )

        if (annotation_gffs_ch != null) {
            MERGE_IR_ISOFORMS(
                annotation_gffs_ch,
                sequence_registry_ch,
                final_decisions_ch,
                merger_script_ch
            )
        }
    }
}
