process FINALIZE_DEEPTMHMM_RESULTS {
    tag 'finalize_deeptmhmm_results'
    publishDir "${params.outdir}/deeptmhmm", mode: 'copy'

    input:
    path initial_manifest
    path generated_results

    output:
    path 'deeptmhmm_results.complete.tsv', emit: manifest, optional: true
    path 'deeptmhmm_results.complete', emit: results_dir, optional: true
    path 'deeptmhmm_results.complete/*', emit: result_dirs, optional: true
    path 'unresolved_deeptmhmm.faa', emit: unresolved_fasta
    path 'deeptmhmm_reuse_summary.tsv', emit: summary

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

    python ${projectDir}/bin/sync_deeptmhmm_dump.py \\
        deeptmhmm_results \\
        ${projectDir}/deeptmhmm_dump

    if grep -q '^>' unresolved_deeptmhmm.faa; then
        echo "WARNING: DeepTMHMM results are incomplete." >&2
        echo "Unresolved sequences were written to:" >&2
        echo "  ${params.outdir}/deeptmhmm/unresolved_deeptmhmm.faa" >&2
        echo "Add completed results under --deeptmhmm_salvage_paths and rerun with -resume." >&2
    else
        mv deeptmhmm_results.tsv deeptmhmm_results.complete.tsv
        mv deeptmhmm_results deeptmhmm_results.complete
    fi
    """
}
