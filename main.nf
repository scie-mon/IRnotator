nextflow.enable.dsl = 2

include { NORMALIZE_PROTEINS } from './modules/local/normalize_input'
include { NORMALIZE_MULTI_GENOME_GFF } from './modules/local/normalize_input'
include { HMMSEARCH_IR } from './modules/local/hmmsearch_ir'
include { COLLECT_HMM_HITS } from './modules/local/collect_hmm_hits'
include { EXTRACT_FASTA_BY_ID } from './modules/local/extract_fasta_by_id'
include { DEEPTMHMM_TOPOLOGY } from './modules/local/deeptmhmm_topology'
include { INITIALIZE_DEEPTMHMM_RESULTS } from './modules/local/initialize_deeptmhmm_results'
include { FINALIZE_DEEPTMHMM_RESULTS } from './modules/local/finalize_deeptmhmm_results'
include { ASSERT_COMPLETE_DEEPTMHMM_RESULTS } from './modules/local/finalize_deeptmhmm_results'
include { COLLECT_DEEPTMHMM_SUMMARY } from './modules/local/collect_deeptmhmm_summary'
include { PARSE_DEEPTMHMM_FEATURES } from './modules/local/parse_deeptmhmm_features'
include { HARDCODED_CLASSIFIER } from './modules/local/hardcoded_classifier'
include { SEED_DECISIONS } from './modules/local/seed_decisions'
include { PREPARE_MANUAL_REVIEW } from './modules/local/prepare_manual_review'
include { GATE_REVIEW_TSV } from './modules/local/gate_review_tsv'
include { APPLY_MANUAL_DECISIONS } from './modules/local/apply_manual_decisions'
include { EXTRACT_PASSED_FASTA } from './modules/local/extract_passed_fasta'
include { MERGE_IR_ISOFORMS as MERGE_UNFILTERED } from './modules/local/merge_ir_isoforms'
include { MERGE_IR_ISOFORMS as MERGE_FILTERED } from './modules/local/merge_ir_isoforms'

params.proteins_faa = params.proteins_faa ?: null
params.genome_fasta = params.genome_fasta ?: null
params.annot_gff = params.annot_gff ?: null
params.gff_protein_attribute = params.gff_protein_attribute ?: 'protein_id'
params.translation_table = params.translation_table ?: 1
params.allow_internal_stops = params.allow_internal_stops ?: false
params.isoform_overlap_fraction = params.isoform_overlap_fraction ?: 0.60
params.hmm_dir = params.hmm_dir ?: 'test/hmms'
params.outdir = params.outdir ?: 'results'
params.deeptmhmm_dir = params.deeptmhmm_dir ?: null
params.run_deeptmhmm = params.run_deeptmhmm == null ? true : params.run_deeptmhmm
params.deeptmhmm_salvage_paths = params.deeptmhmm_salvage_paths ?: []
params.review_tsv = params.review_tsv ?: null

workflow {
    def has_proteins = params.proteins_faa as boolean
    def has_genome = params.genome_fasta as boolean
    def has_gff = params.annot_gff as boolean

    if (has_proteins && has_gff) throw new IllegalArgumentException('The protein+GFF mode is unsupported; use proteins alone or genome+GFF.')
    if (has_proteins && has_genome) throw new IllegalArgumentException('Invalid input: do not provide proteins_faa with genome_fasta.')
    if (!has_proteins && !has_genome) throw new IllegalArgumentException('Invalid input: provide proteins_faa or genome_fasta plus annot_gff.')
    if (has_genome && !has_gff) throw new IllegalArgumentException('Genome FASTA input requires one or more annot_gff files.')
    if (params.isoform_overlap_fraction < 0.0 || params.isoform_overlap_fraction > 1.0) throw new IllegalArgumentException('isoform_overlap_fraction must be between 0 and 1.')

    def review_file = params.review_tsv ? file(params.review_tsv as String) : file("${params.outdir}/manual_review/review.tsv")
    if (params.review_tsv && !review_file.exists()) throw new IllegalArgumentException("Manual review TSV not found: ${review_file.toAbsolutePath()}")

    def review_state = null
    if (has_genome) {
        def merged_dir = file("${params.outdir}/merged_gff")
        def state_file = file("${params.outdir}/merged_gff/review-state.json")
        def command = ['python3', "${projectDir}/bin/resolve_review_state.py", '--merged-dir', merged_dir.toAbsolutePath().toString(), '--output', state_file.toAbsolutePath().toString()]
        if (review_file.exists()) command.addAll(['--review-tsv', review_file.toAbsolutePath().toString()])
        def resolver = command.execute()
        def stdout = new StringBuffer()
        def stderr = new StringBuffer()
        resolver.consumeProcessOutput(stdout, stderr)
        if (resolver.waitFor() != 0) throw new IllegalStateException("Review-state resolution failed: ${stderr}")
        review_state = new groovy.json.JsonSlurper().parse(state_file)
        log.info("Review state: mode=${review_state.mode}; context=${review_state.reviewer_gff}; next=${review_state.target_filtered_gff}")
    }

    normalizer_script_ch = Channel.value(file("${projectDir}/bin/normalize_input.py"))
    multi_normalizer_script_ch = Channel.value(file("${projectDir}/bin/multi_normalize.py"))
    merger_script_ch = Channel.value(file("${projectDir}/bin/merge_ir_isoforms.py"))
    renderer_script_ch = Channel.value(file("${projectDir}/bin/render_gene_cds.py"))

    annotation_gffs_ch = null
    if (has_proteins && !has_genome) {
        normalize_res = NORMALIZE_PROTEINS(Channel.fromPath(params.proteins_faa, checkIfExists: true), normalizer_script_ch)
    } else {
        annotation_gffs_ch = Channel.fromPath(params.annot_gff, checkIfExists: true).collect()
        normalize_res = NORMALIZE_MULTI_GENOME_GFF(Channel.fromPath(params.genome_fasta, checkIfExists: true).collect(), annotation_gffs_ch, normalizer_script_ch, multi_normalizer_script_ch)
    }

    normalized_proteins_ch = normalize_res.normalized.map { proteins_faa, registry_tsv -> proteins_faa }
    sequence_registry_ch = normalize_res.normalized.map { proteins_faa, registry_tsv -> registry_tsv }
    hmm_res = HMMSEARCH_IR(Channel.fromPath("${params.hmm_dir}/*.hmm", checkIfExists: true).combine(normalized_proteins_ch))
    hmm_ids_ch = COLLECT_HMM_HITS(hmm_res.tbl.collect())
    hmm_hit_faa_ch = EXTRACT_FASTA_BY_ID(hmm_ids_ch, normalized_proteins_ch)
    def salvage_paths = (params.deeptmhmm_salvage_paths instanceof Collection ? params.deeptmhmm_salvage_paths.collect { it.toString() } : params.deeptmhmm_salvage_paths.toString().split(',').collect { it.trim() }.findAll { it })
    initial_deeptmhmm_res = INITIALIZE_DEEPTMHMM_RESULTS(hmm_hit_faa_ch, Channel.value(salvage_paths))

    def generated_results_ch
    if (params.run_deeptmhmm) {
        if (!params.deeptmhmm_dir) throw new IllegalArgumentException('Missing deeptmhmm_dir while run_deeptmhmm is enabled.')
        deeptmhmm_res = DEEPTMHMM_TOPOLOGY(initial_deeptmhmm_res.unresolved_fasta.flatten(), Channel.value(file(params.deeptmhmm_dir)))
        generated_results_ch = deeptmhmm_res.results.map { sequence_id, result_dir -> result_dir }.collect().ifEmpty { [] }
    } else {
        generated_results_ch = Channel.value([])
    }

    final_deeptmhmm_res = FINALIZE_DEEPTMHMM_RESULTS(initial_deeptmhmm_res.manifest, generated_results_ch)
    ASSERT_COMPLETE_DEEPTMHMM_RESULTS(final_deeptmhmm_res.unresolved_fasta)
    final_deeptmhmm_dirs_ch = final_deeptmhmm_res.result_dirs.flatten()
    summary_res = COLLECT_DEEPTMHMM_SUMMARY(final_deeptmhmm_res.manifest, final_deeptmhmm_res.results_dir)
    parsed_res = PARSE_DEEPTMHMM_FEATURES(summary_res.summary_dir)
    hardcoded_res = HARDCODED_CLASSIFIER(parsed_res.features)
    seed_res = SEED_DECISIONS(parsed_res.features, hardcoded_res.scores, false)

    def effective_decisions_ch
    def applied_review_tsv_ch
    if (review_file.exists()) {
        def prior_review_dir = file("${params.outdir}/manual_review")
        if (!prior_review_dir.exists()) throw new IllegalArgumentException("Review TSV exists but prior review package is missing: ${prior_review_dir}")
        gate_res = GATE_REVIEW_TSV(review_file.toAbsolutePath().toString(), prior_review_dir)
        effective_decisions_ch = APPLY_MANUAL_DECISIONS(seed_res.decisions, gate_res.review_tsv).decisions
        applied_review_tsv_ch = gate_res.review_tsv
    } else {
        effective_decisions_ch = seed_res.decisions

        def no_review_sentinel = file("${params.outdir}/.no_applied_review.tsv")
        if (!no_review_sentinel.exists()) {
            no_review_sentinel.getParent().toFile().mkdirs()
            no_review_sentinel.toFile().text = ''
        }

        applied_review_tsv_ch = Channel.value(no_review_sentinel)
    }

    EXTRACT_PASSED_FASTA(effective_decisions_ch, hmm_hit_faa_ch)

    if (has_genome) {
        def reviewer_gff_ch
        if (review_state.mode == 'bootstrap') {
            unfiltered_res = MERGE_UNFILTERED(annotation_gffs_ch, sequence_registry_ch, effective_decisions_ch, merger_script_ch, Channel.value('unfiltered'), Channel.value('merged.unfiltered.gff3'), Channel.value('merged.unfiltered.audit.tsv'))
            reviewer_gff_ch = unfiltered_res.merged_gff
        } else {
            reviewer_gff_ch = Channel.fromPath(review_state.reviewer_gff.toString(), checkIfExists: true)
        }

        if (review_state.mode in ['bootstrap', 'advance']) {
            def filtered_name = file(review_state.target_filtered_gff.toString()).getFileName().toString()
            def audit_name = filtered_name.replace('.gff3', '.audit.tsv')
            MERGE_FILTERED(annotation_gffs_ch, sequence_registry_ch, effective_decisions_ch, merger_script_ch, Channel.value('filtered'), Channel.value(filtered_name), Channel.value(audit_name))
        }

        def reviewer_revision = review_state.mode == 'bootstrap' ? -1 : review_state.existing_max_revision
        def written_revision = (review_state.mode == 'bootstrap' ? 0 : (review_state.mode == 'advance' ? review_state.existing_max_revision + 1 : -1))
        PREPARE_MANUAL_REVIEW(
            final_deeptmhmm_dirs_ch.collect(),
            hardcoded_res.scores,
            sequence_registry_ch,
            parsed_res.features,
            effective_decisions_ch,
            applied_review_tsv_ch,
            reviewer_gff_ch,
            file("${projectDir}/assets/reviewer/topology-reviewer.html"),
            renderer_script_ch,
            Channel.value(reviewer_revision),
            Channel.value(written_revision),
        )
    } else {
        log.info('Protein-only mode completed; genome-context manual review is unavailable.')
    }
}
