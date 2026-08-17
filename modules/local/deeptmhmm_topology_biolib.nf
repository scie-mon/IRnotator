process DEEPTMHMM_TOPOLOGY_BIOLIB {
    tag "${proteins_faa.baseName}"
    maxForks 1

    input:
    path proteins_faa
    path runner

    output:
    tuple val(proteins_faa.baseName), path("${proteins_faa.baseName}/"), emit: results

    script:
    """
    mkdir -p "${proteins_faa.baseName}"

    set +e
    python "${runner}" \
        "${proteins_faa}" \
        --outdir biolib_out \
        --app "${params.deeptmhmm_model}" \
        --manifest biolib_manifest.tsv
    set -e

    if [[ -d "biolib_out/${proteins_faa.baseName}" ]] && compgen -G "biolib_out/${proteins_faa.baseName}/*" > /dev/null; then
        mv "biolib_out/${proteins_faa.baseName}"/* "${proteins_faa.baseName}/"
    fi
    """
}