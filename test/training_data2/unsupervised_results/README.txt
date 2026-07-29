Unsupervised topology analysis outputs
====================================

Accepted input:
- a directory of one-sequence .3line files
- or one multi-sequence .3line file

Main files:
- dataset_summary.csv: overall dataset stats after parsing and cropping
- sequence_assignments.csv: one row per sequence with crop results, UMAP coordinates, and cluster IDs
- clustering_metrics.csv: silhouette scores for clustering runs
- cluster_summaries.csv: cluster sizes, mean lengths, and failed-overlay enrichment
- cluster_top_features.csv: most enriched k-mer features per cluster
- umap_*.png: 2D UMAP visualizations for each clustering run
- *_failed_overlay.png: same UMAP plots with failed PNG-linked sequences marked by x
- dendrogram_subsample.png: hierarchical clustering dendrogram on all or a reproducible subsample
- length_distributions.png, raw_vs_cropped_lengths_scatter.png: crop-effect diagnostics

Crop rule:
- remove leading S
- then remove leading O
- then remove trailing I until the last character is not I
