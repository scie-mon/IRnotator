process DEEPTMHMM_TOPOLOGY {
    tag "${proteins_faa.baseName}"

    publishDir "${params.outdir}/deeptmhmm", mode: 'copy'

    input:
    path proteins_faa

    output:
    path "biolib_results", emit: results

    script:
    def cmd
    if( params.deeptmhmm_mode == 'local' ) {
        cmd = "biolib run --local '${params.deeptmhmm_model}' --fasta ${proteins_faa}"
    }
    else if( params.deeptmhmm_mode == 'local_gpu' ) {
        cmd = "export BIOLIB_DOCKER_RUNTIME=nvidia && biolib run --local '${params.deeptmhmm_model}' --fasta ${proteins_faa}"
    }
    else if( params.deeptmhmm_mode == 'cloud' ) {
        cmd = "biolib run '${params.deeptmhmm_model}' --fasta ${proteins_faa}"
    }
    else {
        throw new IllegalArgumentException("Unsupported --deeptmhmm_mode: ${params.deeptmhmm_mode}")
    }

    """
    ${cmd}
    """
}
