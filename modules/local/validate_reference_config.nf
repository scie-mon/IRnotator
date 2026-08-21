process VALIDATE_REFERENCE_CONFIG {
    tag 'reference configuration'

    input:
    val hmm_dir
    val salvage_paths

    output:
    stdout emit: audit

    script:
    def paths = salvage_paths instanceof Collection ? salvage_paths : [salvage_paths]
    def salvage_args = paths
        .findAll { value -> value != null && value.toString().trim() }
        .collect { value -> "--salvage-path '${value.toString().replace("'", "'\\\"'\\\"'")}'" }
        .join(' ')
    def project_dir = projectDir.toString().replace("'", "'\\\"'\\\"'")
    def validator = "${projectDir}/bin/validate_reference_config.py".replace("'", "'\\\"'\\\"'")
    def configured_hmm_dir = hmm_dir.toString().replace("'", "'\\\"'\\\"'")

    """
    cd '${project_dir}'
    python '${validator}' \\
        --project-dir '${project_dir}' \\
        --hmm-dir '${configured_hmm_dir}' \\
        ${salvage_args}
    """
}
