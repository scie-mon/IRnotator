process PARSE_DEEPTMHMM_FEATURES {
    tag "parse_deeptmhmm"

    input:
    path deeptmhmm_dir

    output:
    path "deeptmhmm_features.tsv", emit: features

    script:
    """
    python ${projectDir}/bin/parse_deeptmhmm_features.py \
        ${deeptmhmm_dir}/TMRs.gff3 \
        deeptmhmm_features.tsv
    """
}