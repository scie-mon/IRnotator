Supervised DeepTMHMM topology CNN
=================================

Labels are inferred from the PNG sister directory:
- PNGs in root are labeled pass
- PNGs in failed/ are labeled false_positive
- .3line files are matched by shared basename

Training outputs:
- model.pt: trained PyTorch model weights
- model_metadata.json: preprocessing and architecture parameters
- training_predictions_oof.csv: out-of-fold predictions for all labeled examples
- cv_fold_metrics.csv: average precision per fold
- threshold_scan.csv: precision/recall across thresholds
- training_summary.csv: overall summary and confusion matrix

Prediction rule:
- score >= decision_threshold => false_positive
- else => pass
