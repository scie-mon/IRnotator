process ECHO_INPUT {
    tag "${input_file.baseName}"

    publishDir "${params.outdir}/echo", mode: 'copy'

    input:
    path input_file

    output:
    path "echoed_${input_file.getName()}"

    script:
    """
    cp ${input_file} echoed_${input_file.getName()}
    echo "Processed: ${input_file.getName()}" >> echoed_${input_file.getName()}
    """
}
