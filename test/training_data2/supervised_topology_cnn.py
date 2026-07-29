#!/usr/bin/env python3
import argparse
import json
import random
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score, precision_recall_curve, confusion_matrix


ALPHABET = ['P', 'S', 'O', 'M', 'I', 'X']
CHAR_TO_IDX = {c: i for i, c in enumerate(ALPHABET)}
PAD_CHAR = 'X'
PAD_IDX = CHAR_TO_IDX[PAD_CHAR]


def parse_args():
    p = argparse.ArgumentParser(
        description="Train or apply a supervised CNN on DeepTMHMM topology strings using PNG failed/ as labels."
    )
    sub = p.add_subparsers(dest="command", required=True)

    tr = sub.add_parser("train", help="Train model from .3line directory + PNG sister dir labels")
    tr.add_argument("input_3line_dir", help="Directory containing one-sequence .3line files")
    tr.add_argument("--png-sister-dir", required=True,
                    help="Directory containing PNGs with same basenames and failed/ subdir for negative labels")
    tr.add_argument("--outdir", default="supervised_topology_cnn_results")
    tr.add_argument("--epochs", type=int, default=30)
    tr.add_argument("--batch-size", type=int, default=64)
    tr.add_argument("--lr", type=float, default=1e-3)
    tr.add_argument("--embed-dim", type=int, default=16)
    tr.add_argument("--channels", type=int, default=64)
    tr.add_argument("--kernel-sizes", default="3,5,7")
    tr.add_argument("--dropout", type=float, default=0.25)
    tr.add_argument("--n-splits", type=int, default=5)
    tr.add_argument("--seed", type=int, default=42)
    tr.add_argument("--decision-threshold", type=float, default=None,
                    help="Optional fixed threshold for fail probability; if omitted, choose CV threshold targeting high precision")
    tr.add_argument("--min-positive-png", type=int, default=1,
                    help="Minimum number of non-failed PNGs required to include a sample as positive label")

    pr = sub.add_parser("predict", help="Apply trained model to .3line directory")
    pr.add_argument("input_3line_dir", help="Directory containing one-sequence .3line files")
    pr.add_argument("--model-dir", required=True, help="Directory produced by train command")
    pr.add_argument("--outdir", default="supervised_topology_predictions")

    return p.parse_args()


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def crop_topology(s: str) -> str:
    s = str(s).strip().upper()
    i = 0
    while i < len(s) and s[i] == 'S':
        i += 1
    while i < len(s) and s[i] == 'O':
        i += 1
    j = len(s)
    while j > i and s[j - 1] == 'I':
        j -= 1
    return s[i:j]


def parse_3line_file(path: Path):
    lines = [line.rstrip("\n") for line in path.read_text().splitlines() if line.strip()]
    if len(lines) % 3 != 0:
        raise ValueError(f"Input file does not contain a multiple of 3 non-empty lines: {path}")
    if len(lines) != 3:
        raise ValueError(f"Expected one sequence per .3line file, got {len(lines)//3} records in {path}")
    header, aa, topo = lines
    if not header.startswith(">"):
        raise ValueError(f"Expected header line starting with '>' in {path}")
    header_text = header[1:].strip()
    if " | " in header_text:
        header_seq_id, prediction_label = header_text.split(" | ", 1)
    else:
        header_seq_id, prediction_label = header_text, ""
    return {
        "id": path.stem,
        "header_seq_id": header_seq_id,
        "prediction_label": prediction_label,
        "header_text": header_text,
        "aa_sequence": aa.strip(),
        "topology_string": topo.strip().upper(),
        "source_file": str(path),
    }


def load_3line_dir(input_dir: Path) -> pd.DataFrame:
    files = sorted(list(input_dir.glob("*.3line")) + list(input_dir.glob("*.txt")) + list(input_dir.glob("*.3line.txt")))
    if not files:
        raise ValueError(f"No .3line-like files found in directory: {input_dir}")
    records = [parse_3line_file(f) for f in files]
    df = pd.DataFrame(records)
    if df["id"].duplicated().any():
        dupes = df.loc[df["id"].duplicated(), "id"].tolist()[:10]
        raise ValueError(f"Duplicate file-derived IDs detected. Example duplicates: {dupes}")
    df["cropped_topology"] = df["topology_string"].apply(crop_topology)
    df["raw_len"] = df["topology_string"].str.len()
    df["cropped_len"] = df["cropped_topology"].str.len()
    return df


def collect_png_labels(png_sister_dir: Path, min_positive_png: int = 1) -> pd.DataFrame:
    failed_dir = png_sister_dir / "failed"
    all_png = [p for p in png_sister_dir.glob("*.png") if p.is_file()]
    failed_png = list(failed_dir.glob("*.png")) if failed_dir.exists() else []

    pos_ids = sorted({p.stem for p in all_png} - {p.stem for p in failed_png})
    neg_ids = sorted({p.stem for p in failed_png})

    rows = []
    for sid in pos_ids:
        rows.append({"id": sid, "label": 0, "label_name": "pass", "png_label_source": "png_root"})
    for sid in neg_ids:
        rows.append({"id": sid, "label": 1, "label_name": "false_positive", "png_label_source": "failed_subdir"})

    df = pd.DataFrame(rows).drop_duplicates(subset=["id"], keep="last")
    if len(df) == 0:
        raise ValueError(f"No PNG-derived labels found in {png_sister_dir}")
    return df


class TopologyDataset(Dataset):
    def __init__(self, encoded: np.ndarray, labels: np.ndarray | None = None):
        self.encoded = torch.tensor(encoded, dtype=torch.long)
        self.labels = None if labels is None else torch.tensor(labels, dtype=torch.float32)

    def __len__(self):
        return len(self.encoded)

    def __getitem__(self, idx):
        if self.labels is None:
            return self.encoded[idx]
        return self.encoded[idx], self.labels[idx]


class TopologyCNN(nn.Module):
    def __init__(self, vocab_size: int, embed_dim: int, channels: int, kernel_sizes: List[int], dropout: float):
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=PAD_IDX)
        self.convs = nn.ModuleList([
            nn.Conv1d(embed_dim, channels, kernel_size=k, padding=k // 2)
            for k in kernel_sizes
        ])
        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(channels * len(kernel_sizes), 1)

    def forward(self, x):
        x = self.embedding(x)
        x = x.transpose(1, 2)
        pooled = []
        for conv in self.convs:
            y = torch.relu(conv(x))
            y = torch.max(y, dim=2).values
            pooled.append(y)
        z = torch.cat(pooled, dim=1)
        z = self.dropout(z)
        return self.fc(z).squeeze(1)


def encode_topologies(topologies: List[str], max_len: int) -> np.ndarray:
    arr = np.full((len(topologies), max_len), PAD_IDX, dtype=np.int64)
    for i, s in enumerate(topologies):
        s = ''.join(ch if ch in CHAR_TO_IDX else PAD_CHAR for ch in s)
        s = s[:max_len]
        idxs = [CHAR_TO_IDX[ch] for ch in s]
        arr[i, :len(idxs)] = idxs
    return arr


def train_one_fold(model, train_loader, val_loader, device, epochs, lr, pos_weight):
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([pos_weight], dtype=torch.float32, device=device))
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    best_state = None
    best_ap = -1.0
    history = []

    for epoch in range(1, epochs + 1):
        model.train()
        train_losses = []
        for xb, yb in train_loader:
            xb = xb.to(device)
            yb = yb.to(device)
            optimizer.zero_grad()
            logits = model(xb)
            loss = criterion(logits, yb)
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())

        model.eval()
        val_true, val_score = [], []
        with torch.no_grad():
            for xb, yb in val_loader:
                xb = xb.to(device)
                logits = model(xb)
                probs = torch.sigmoid(logits).cpu().numpy()
                val_score.extend(probs.tolist())
                val_true.extend(yb.numpy().tolist())

        val_true = np.asarray(val_true, dtype=int)
        val_score = np.asarray(val_score, dtype=float)
        val_ap = average_precision_score(val_true, val_score)
        row = {"epoch": epoch, "train_loss": float(np.mean(train_losses)), "val_ap": float(val_ap)}
        history.append(row)
        if val_ap > best_ap:
            best_ap = val_ap
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    model.load_state_dict(best_state)
    return model, pd.DataFrame(history), best_ap


def choose_threshold(y_true: np.ndarray, y_score: np.ndarray) -> Tuple[float, pd.DataFrame]:
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    rows = []
    best_thr = 0.5
    best_recall = -1.0
    for i, thr in enumerate(thresholds):
        p = float(precision[i + 1])
        r = float(recall[i + 1])
        rows.append({"threshold": float(thr), "precision": p, "recall": r})
        if p >= 0.95 and r > best_recall:
            best_recall = r
            best_thr = float(thr)
    if best_recall < 0:
        f1_best = -1.0
        for row in rows:
            p, r = row["precision"], row["recall"]
            f1 = 0.0 if (p + r) == 0 else (2 * p * r / (p + r))
            if f1 > f1_best:
                f1_best = f1
                best_thr = row["threshold"]
    return best_thr, pd.DataFrame(rows)


def predict_scores(model, loader, device):
    model.eval()
    scores = []
    with torch.no_grad():
        for batch in loader:
            xb = batch[0] if isinstance(batch, (tuple, list)) else batch
            xb = xb.to(device)
            logits = model(xb)
            probs = torch.sigmoid(logits).cpu().numpy()
            scores.extend(probs.tolist())
    return np.asarray(scores, dtype=float)


def run_train(args):
    set_seed(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    df = load_3line_dir(Path(args.input_3line_dir))
    label_df = collect_png_labels(Path(args.png_sister_dir), min_positive_png=args.min_positive_png)
    train_df = df.merge(label_df, on="id", how="inner")

    if len(train_df) == 0:
        raise ValueError("No overlap between .3line basenames and PNG-derived labels")
    if train_df["label"].nunique() < 2:
        raise ValueError("Need both pass and false_positive labels for training")

    max_len = int(train_df["cropped_len"].max())
    X = encode_topologies(train_df["cropped_topology"].tolist(), max_len=max_len)
    y = train_df["label"].to_numpy(dtype=int)

    kernel_sizes = [int(x) for x in args.kernel_sizes.split(",") if x.strip()]
    skf = StratifiedKFold(n_splits=args.n_splits, shuffle=True, random_state=args.seed)

    all_oof_scores = np.zeros(len(train_df), dtype=float)
    fold_metrics = []
    history_frames = []

    for fold, (tr_idx, va_idx) in enumerate(skf.split(X, y), start=1):
        X_tr, X_va = X[tr_idx], X[va_idx]
        y_tr, y_va = y[tr_idx], y[va_idx]

        train_loader = DataLoader(TopologyDataset(X_tr, y_tr), batch_size=args.batch_size, shuffle=True)
        val_loader = DataLoader(TopologyDataset(X_va, y_va), batch_size=args.batch_size, shuffle=False)

        n_pos = int(y_tr.sum())
        n_neg = int((y_tr == 0).sum())
        pos_weight = n_neg / max(n_pos, 1)

        model = TopologyCNN(
            vocab_size=len(ALPHABET),
            embed_dim=args.embed_dim,
            channels=args.channels,
            kernel_sizes=kernel_sizes,
            dropout=args.dropout,
        ).to(device)

        model, hist_df, best_ap = train_one_fold(model, train_loader, val_loader, device, args.epochs, args.lr, pos_weight)
        hist_df.insert(0, "fold", fold)
        history_frames.append(hist_df)

        val_scores = predict_scores(model, val_loader, device)
        all_oof_scores[va_idx] = val_scores
        fold_ap = average_precision_score(y_va, val_scores)
        fold_metrics.append({"fold": fold, "n_train": len(tr_idx), "n_val": len(va_idx), "val_ap": float(fold_ap)})

    if args.decision_threshold is None:
        threshold, threshold_df = choose_threshold(y, all_oof_scores)
    else:
        threshold = float(args.decision_threshold)
        precision, recall, thresholds = precision_recall_curve(y, all_oof_scores)
        threshold_df = pd.DataFrame({
            "threshold": thresholds,
            "precision": precision[1:],
            "recall": recall[1:],
        })

    oof_pred = (all_oof_scores >= threshold).astype(int)
    cm = confusion_matrix(y, oof_pred, labels=[0, 1])

    train_df = train_df.copy()
    train_df["oof_false_positive_score"] = all_oof_scores
    train_df["oof_pred_label"] = oof_pred
    train_df["oof_pred_name"] = np.where(oof_pred == 1, "false_positive", "pass")

    full_loader = DataLoader(TopologyDataset(X, y), batch_size=args.batch_size, shuffle=True)
    n_pos_full = int(y.sum())
    n_neg_full = int((y == 0).sum())
    pos_weight_full = n_neg_full / max(n_pos_full, 1)

    final_model = TopologyCNN(
        vocab_size=len(ALPHABET),
        embed_dim=args.embed_dim,
        channels=args.channels,
        kernel_sizes=kernel_sizes,
        dropout=args.dropout,
    ).to(device)
    final_model, final_hist, final_ap = train_one_fold(final_model, full_loader, full_loader, device, args.epochs, args.lr, pos_weight_full)

    torch.save(final_model.state_dict(), outdir / "model.pt")

    metadata = {
        "alphabet": ALPHABET,
        "pad_char": PAD_CHAR,
        "pad_idx": PAD_IDX,
        "max_len": max_len,
        "embed_dim": args.embed_dim,
        "channels": args.channels,
        "kernel_sizes": kernel_sizes,
        "dropout": args.dropout,
        "decision_threshold": threshold,
        "label_definition": {"0": "pass", "1": "false_positive"},
        "device_used_for_training": str(device),
        "cropping_rule": "remove leading S, then leading O, then trailing I until last char is not I",
    }
    (outdir / "model_metadata.json").write_text(json.dumps(metadata, indent=2))

    pd.DataFrame(fold_metrics).to_csv(outdir / "cv_fold_metrics.csv", index=False)
    pd.concat(history_frames, ignore_index=True).to_csv(outdir / "cv_training_history.csv", index=False)
    threshold_df.to_csv(outdir / "threshold_scan.csv", index=False)
    train_df.to_csv(outdir / "training_predictions_oof.csv", index=False)
    final_hist.to_csv(outdir / "final_model_training_history.csv", index=False)

    summary = pd.DataFrame([{
        "n_sequences_used": int(len(train_df)),
        "n_pass": int((train_df['label'] == 0).sum()),
        "n_false_positive": int((train_df['label'] == 1).sum()),
        "max_cropped_len": int(max_len),
        "cv_average_precision": float(average_precision_score(y, all_oof_scores)),
        "chosen_threshold": float(threshold),
        "tn": int(cm[0, 0]),
        "fp": int(cm[0, 1]),
        "fn": int(cm[1, 0]),
        "tp": int(cm[1, 1]),
        "device_used_for_training": str(device),
    }])
    summary.to_csv(outdir / "training_summary.csv", index=False)

    readme = outdir / "README.txt"
    readme.write_text(
        "Supervised DeepTMHMM topology CNN\n"
        "=================================\n\n"
        "Labels are inferred from the PNG sister directory:\n"
        "- PNGs in root are labeled pass\n"
        "- PNGs in failed/ are labeled false_positive\n"
        "- .3line files are matched by shared basename\n\n"
        "Training outputs:\n"
        "- model.pt: trained PyTorch model weights\n"
        "- model_metadata.json: preprocessing and architecture parameters\n"
        "- training_predictions_oof.csv: out-of-fold predictions for all labeled examples\n"
        "- cv_fold_metrics.csv: average precision per fold\n"
        "- threshold_scan.csv: precision/recall across thresholds\n"
        "- training_summary.csv: overall summary and confusion matrix\n\n"
        "Prediction rule:\n"
        "- score >= decision_threshold => false_positive\n"
        "- else => pass\n"
    )



def load_model(model_dir: Path, device):
    metadata = json.loads((model_dir / "model_metadata.json").read_text())
    model = TopologyCNN(
        vocab_size=len(metadata["alphabet"]),
        embed_dim=metadata["embed_dim"],
        channels=metadata["channels"],
        kernel_sizes=metadata["kernel_sizes"],
        dropout=metadata["dropout"],
    ).to(device)
    state = torch.load(model_dir / "model.pt", map_location=device)
    model.load_state_dict(state)
    model.eval()
    return model, metadata


def run_predict(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    df = load_3line_dir(Path(args.input_3line_dir))
    model, metadata = load_model(Path(args.model_dir), device)
    X = encode_topologies(df["cropped_topology"].tolist(), max_len=int(metadata["max_len"]))
    loader = DataLoader(TopologyDataset(X, None), batch_size=128, shuffle=False)
    scores = predict_scores(model, loader, device)

    threshold = float(metadata["decision_threshold"])
    pred = (scores >= threshold).astype(int)

    df = df.copy()
    df["false_positive_score"] = scores
    df["pred_label"] = pred
    df["pred_name"] = np.where(pred == 1, "false_positive", "pass")
    df["decision_threshold"] = threshold
    df.to_csv(outdir / "predictions.csv", index=False)

    summary = pd.DataFrame([{
        "n_sequences": int(len(df)),
        "n_pred_pass": int((pred == 0).sum()),
        "n_pred_false_positive": int((pred == 1).sum()),
        "decision_threshold": threshold,
        "device_used_for_prediction": str(device),
    }])
    summary.to_csv(outdir / "prediction_summary.csv", index=False)

    (outdir / "README.txt").write_text(
        "Prediction outputs\n"
        "==================\n\n"
        "- predictions.csv: per-sequence scores and calls\n"
        "- prediction_summary.csv: counts by predicted label\n"
    )


if __name__ == "__main__":
    args = parse_args()
    if args.command == "train":
        run_train(args)
    elif args.command == "predict":
        run_predict(args)
