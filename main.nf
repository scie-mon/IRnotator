nextflow.enable.dsl=2

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
include { EXTRACT_PASSED_FASTA } from './modules/local/extract_passed_fasta'

params.proteins_faa  = params.proteins_faa  ?: 'test/proteins.faa'
params.hmm_dir       = params.hmm_dir       ?: 'test/hmms'
params.outdir        = params.outdir        ?: 'results'
params.deeptmhmm_dir = params.deeptmhmm_dir ?: null
params.manual_review = params.manual_review ?: true
params.reviewer_html = params.reviewer_html ?: "${projectDir}/assets/reviewer/topology-reviewer.html"

workflow {
    proteins_ch      = Channel.fromPath(params.proteins_faa, checkIfExists: true)
    hmms_ch          = Channel.fromPath("${params.hmm_dir}/*.hmm", checkIfExists: true)
    deeptmhmm_dir_ch = Channel.value(file(params.deeptmhmm_dir))

    hmm_input_ch = hmms_ch.combine(proteins_ch)

    hmm_res     = HMMSEARCH_IR(hmm_input_ch)
    hmm_tbls_ch = hmm_res.tbl

    hmm_ids_ch     = COLLECT_HMM_HITS(hmm_tbls_ch.collect())
    hmm_hit_faa_ch = EXTRACT_FASTA_BY_ID(hmm_ids_ch, proteins_ch)

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

    if ( params.manual_review ) {
        prepare_res = PREPARE_MANUAL_REVIEW(
            deeptmhmm_res.results.map { id, dir -> dir }.collect(),
            hardcoded_res.scores,
            file(params.reviewer_html)
        )

        def review_tsv_path = file("${params.outdir}/manual_review/review.tsv").toAbsolutePath().toString()

        gate_res = GATE_REVIEW_TSV(
            review_tsv_path,
            prepare_res.review_dir
        )

        decisions_ch = APPLY_MANUAL_DECISIONS(
            seed_res.decisions,
            gate_res.review_tsv
        ).decisions
        passed_res = EXTRACT_PASSED_FASTA(
            decisions_ch,
            hmm_hit_faa_ch
        )
        // passed_res.fasta

    } else {
        decisions_ch = seed_res.decisions
        passed_res = EXTRACT_PASSED_FASTA(
            decisions_ch,
            hmm_hit_faa_ch
        )
        // passed_res.fasta
    }

    // decisions_ch → next: FASTA emit (final_decision == accept)
}
