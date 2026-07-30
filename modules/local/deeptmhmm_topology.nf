process DEEPTMHMM_TOPOLOGY {
    tag "${proteins_faa.baseName}"
    maxForks 1
    publishDir "${params.outdir}/deeptmhmm", mode: 'copy'

    input:
    path proteins_faa
    path deeptmhmm_dir, stageAs: 'deeptmhmm'

    output:
    tuple val(proteins_faa.baseName), path("${proteins_faa.baseName}/"), emit: results

    script:
    if (params.deeptmhmm_mode == 'container') {
        """
        rm -rf deeptmhmm_tmp ${proteins_faa.baseName}
        cp ${proteins_faa} input.faa

        FASTA_PATH=\$(readlink -f input.faa)
        TMP_OUT=\$(pwd)/deeptmhmm_tmp
        FINAL_OUT=\$(pwd)/${proteins_faa.baseName}

        cd deeptmhmm
        python predict.py \\
            --fasta "\${FASTA_PATH}" \\
            --output-dir "\${TMP_OUT}"

        cd ..
        mkdir -p "\${FINAL_OUT}"
        mv deeptmhmm_tmp/* "\${FINAL_OUT}/"
        rmdir deeptmhmm_tmp
        """
    }
    else {
        throw new IllegalArgumentException(
            "Unsupported --deeptmhmm_mode: ${params.deeptmhmm_mode}. Use 'container' for local per-sequence execution."
        )
    }
}
