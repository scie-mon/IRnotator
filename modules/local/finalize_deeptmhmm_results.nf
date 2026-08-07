process FINALIZE_DEEPTMHMM_RESULTS {
    tag 'finalize_deeptmhmm_results'
    publishDir "${params.outdir}/deeptmhmm", mode: 'copy'

    input:
    path initial_manifest
    path generated_results

    output:
    path 'deeptmhmm_results.tsv', emit: manifest
    path 'deeptmhmm_reuse_summary.tsv', emit: summary
    path 'deeptmhmm_results', emit: results_dir
    path 'deeptmhmm_results/*', emit: result_dirs
    path 'unresolved_deeptmhmm.faa', emit: unresolved_fasta

    script:
    def generated_paths = generated_results.collect { it.toString() }.join('\n')
    """
    cat > generated_results.txt <<'EOF'
${generated_paths}
EOF

    python ${projectDir}/bin/finalize_deeptmhmm_results.py \\
        ${initial_manifest} \\
        generated_results.txt \\
        deeptmhmm_results.tsv \\
        deeptmhmm_results \\
        unresolved_deeptmhmm.faa \\
        deeptmhmm_reuse_summary.tsv
    """
}

process ASSERT_COMPLETE_DEEPTMHMM_RESULTS {
    tag 'assert_complete_deeptmhmm_results'
    input:
    path unresolved_fasta
    script:
    """
    if grep -q '^>' ${unresolved_fasta}; then
        echo 'DeepTMHMM results are missing for one or more sequences:' >&2
        cat ${unresolved_fasta} >&2
        exit 1
    fi
    """
}
