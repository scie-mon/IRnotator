nextflow.enable.dsl=2

include { HMMSEARCH_IR }              from './modules/local/hmmsearch_ir'
include { COLLECT_HMM_HITS }          from './modules/local/collect_hmm_hits'
include { EXTRACT_FASTA_BY_ID }       from './modules/local/extract_fasta_by_id'
include { DEEPTMHMM_TOPOLOGY }        from './modules/local/deeptmhmm_topology'
include { PARSE_DEEPTMHMM_FEATURES }  from './modules/local/parse_deeptmhmm_features'
include { HARDCODED_CLASSIFIER }      from './modules/local/hardcoded_classifier'
include { SEED_DECISIONS }            from './modules/local/seed_decisions'

params.proteins_faa = params.proteins_faa ?: 'test/proteins.faa'
params.hmm_dir      = params.hmm_dir ?: 'test/hmms'
params.outdir       = params.outdir ?: 'results'

workflow {
    proteins_ch = Channel.fromPath(params.proteins_faa, checkIfExists: true)
    hmms_ch     = Channel.fromPath("${params.hmm_dir}/*.hmm", checkIfExists: true)

    hmm_input_ch = hmms_ch.combine(proteins_ch)

    hmm_res = HMMSEARCH_IR(hmm_input_ch)
    hmm_tbls_ch = hmm_res.tbl

    hmm_ids_ch = COLLECT_HMM_HITS(hmm_tbls_ch.collect())
    hmm_hit_faa_ch = EXTRACT_FASTA_BY_ID(hmm_ids_ch, proteins_ch)

    deeptmhmm_res = DEEPTMHMM_TOPOLOGY(hmm_hit_faa_ch)

    parsed_res = PARSE_DEEPTMHMM_FEATURES(deeptmhmm_res.results)
    hardcoded_res = HARDCODED_CLASSIFIER(parsed_res.features)
    SEED_DECISIONS(parsed_res.features, hardcoded_res.scores)
}
