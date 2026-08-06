process HARDCODED_CLASSIFIER {
    tag "hardcoded_classifier"

    input:
    path features_tsv

    output:
    path "hardcoded_scores.tsv", emit: scores

    script:
    """
    python ${projectDir}/bin/hardcoded_classifier.py \
        ${features_tsv} \
        hardcoded_scores.tsv
    """
}
