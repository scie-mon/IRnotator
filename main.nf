nextflow.enable.dsl=2

include { ECHO_INPUT } from './modules/local/echo_input'

params.input  = params.input ?: 'assets/test_input.txt'
params.outdir = params.outdir ?: 'results'

workflow {
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    ECHO_INPUT(input_ch)
}
