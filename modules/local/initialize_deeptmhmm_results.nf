process INITIALIZE_DEEPTMHMM_RESULTS {
    tag 'initialize_deeptmhmm_results'
    cache false

    input:
    path hmm_hit_faa
    val salvage_paths

    output:
    path 'deeptmhmm_results.initial.tsv', emit: manifest
    path 'deeptmhmm_initial_summary.tsv', emit: summary
    path 'unmatched/*.faa', emit: unresolved_fasta, optional: true

    script:
    def roots_text = salvage_paths.join('\n')
    """
    cat > salvage_roots.txt <<'EOF'
${roots_text}
EOF

    python ${projectDir}/bin/initialize_deeptmhmm_results.py \\
        ${hmm_hit_faa} \\
        salvage_roots.txt \\
        deeptmhmm_results.initial.tsv \\
        unmatched \\
        deeptmhmm_initial_summary.tsv
    """
}
