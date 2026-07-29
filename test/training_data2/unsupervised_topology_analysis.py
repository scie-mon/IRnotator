#!/usr/bin/env python3
import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from sklearn.feature_extraction import DictVectorizer
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans, AgglomerativeClustering, DBSCAN
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram
import umap.umap_ as umap


def parse_args():
    p = argparse.ArgumentParser(
        description="Unsupervised analysis of DeepTMHMM topology strings from either a directory of one-sequence .3line files or a single multi-sequence .3line file."
    )
    p.add_argument("input_path", help="Path to a directory of .3line files or one multi-sequence .3line file")
    p.add_argument("--outdir", default="unsupervised_topology_results", help="Output directory")
    p.add_argument("--png-sister-dir", default=None,
                   help="Optional directory containing PNGs with same basenames/IDs, plus optional failed/ subdir")
    p.add_argument("--k-values", default="2,3,4", help="Comma-separated k values for k-means and agglomerative clustering")
    p.add_argument("--umap-neighbors", type=int, default=15)
    p.add_argument("--umap-min-dist", type=float, default=0.1)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--max-dendrogram", type=int, default=1200,
                   help="Maximum number of sequences to include in dendrogram plot")
    return p.parse_args()


def parse_3line_file(path: Path):
    lines = [line.rstrip("\n") for line in path.read_text().splitlines() if line.strip()]
    if len(lines) % 3 != 0:
        raise ValueError(f"Input file does not contain a multiple of 3 non-empty lines: {path}")
    records = []
    multi_record = (len(lines) // 3) > 1
    for i in range(0, len(lines), 3):
        header, aa, topo = lines[i:i+3]
        if not header.startswith(">"):
            raise ValueError(f"Expected header line starting with '>' in {path} at logical line {i+1}")
        header_text = header[1:].strip()
        if " | " in header_text:
            header_seq_id, prediction_label = header_text.split(" | ", 1)
        else:
            header_seq_id, prediction_label = header_text, ""
        if multi_record:
            seq_id = f"{path.stem}__{(i // 3) + 1}"
        else:
            seq_id = path.stem
        records.append({
            "id": seq_id,
            "header_seq_id": header_seq_id,
            "prediction_label": prediction_label,
            "header_text": header_text,
            "aa_sequence": aa.strip(),
            "topology_string": topo.strip().upper(),
            "source_file": str(path)
        })
    return records


def load_input(input_path: Path):
    if input_path.is_dir():
        files = sorted(list(input_path.glob("*.3line")) + list(input_path.glob("*.txt")) + list(input_path.glob("*.3line.txt")))
        if not files:
            raise ValueError(f"No .3line-like files found in directory: {input_path}")
        records = []
        for f in files:
            records.extend(parse_3line_file(f))
        return pd.DataFrame(records), files
    elif input_path.is_file():
        return pd.DataFrame(parse_3line_file(input_path)), [input_path]
    else:
        raise ValueError(f"Input path does not exist: {input_path}")


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


def kmers(seq, k):
    if len(seq) < k:
        return []
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]


def seq_to_features(seq, ks=(2, 3)):
    c = Counter()
    c["cropped_len"] = len(seq)
    for ch in set(seq):
        c[f"char_{ch}"] = seq.count(ch)
    for k in ks:
        for token in kmers(seq, k):
            c[f"k{k}_{token}"] += 1
    return dict(c)


def detect_failed_ids(png_sister_dir: Path | None):
    if png_sister_dir is None:
        return set()
    failed_dir = png_sister_dir / "failed"
    if not failed_dir.exists() or not failed_dir.is_dir():
        return set()
    return {p.stem for p in failed_dir.glob("*.png")}


def summarize_clusters(df, col):
    rows = []
    for cluster_id, sub in df.groupby(col):
        rows.append({
            "cluster_method": col,
            "cluster_id": cluster_id,
            "n_sequences": len(sub),
            "mean_cropped_len": round(sub["cropped_len"].mean(), 3),
            "median_cropped_len": round(sub["cropped_len"].median(), 3),
            "n_failed_overlay": int(sub["is_failed_overlay"].sum()),
            "failed_overlay_fraction": round(float(sub["is_failed_overlay"].mean()), 5),
        })
    return pd.DataFrame(rows).sort_values(["cluster_method", "cluster_id"])


def top_features_per_cluster(X, feature_names, labels, top_n=15):
    rows = []
    labels = np.asarray(labels)
    unique = sorted(pd.unique(labels))
    for cl in unique:
        mask = labels == cl
        if mask.sum() == 0:
            continue
        centroid = X[mask].mean(axis=0)
        top_idx = np.argsort(centroid)[::-1][:top_n]
        for rank, idx in enumerate(top_idx, start=1):
            rows.append({
                "cluster_id": cl,
                "rank": rank,
                "feature": feature_names[idx],
                "mean_value": float(centroid[idx]),
            })
    return pd.DataFrame(rows)


def make_umap_plot(df, x_col, y_col, color_col, title, outpath, failed_overlay=False):
    plt.figure(figsize=(9, 7))
    cats = sorted(pd.unique(df[color_col]))
    cmap = plt.cm.get_cmap("tab10", max(len(cats), 3))
    color_map = {cat: cmap(i) for i, cat in enumerate(cats)}

    if failed_overlay:
        normal = df[~df["is_failed_overlay"]]
        failed = df[df["is_failed_overlay"]]
        plt.scatter(normal[x_col], normal[y_col],
                    c=[color_map[v] for v in normal[color_col]], s=20, alpha=0.75, linewidths=0)
        if len(failed) > 0:
            plt.scatter(failed[x_col], failed[y_col],
                        c=[color_map[v] for v in failed[color_col]],
                        s=80, alpha=0.98, marker='x', linewidths=1.6)
        legend1 = [Line2D([0], [0], marker='o', color='w', label=str(cat),
                          markerfacecolor=color_map[cat], markersize=8) for cat in cats]
        legend2 = [
            Line2D([0], [0], marker='o', color='black', label='normal', linestyle='None', markersize=6),
            Line2D([0], [0], marker='x', color='black', label='failed overlay', linestyle='None', markersize=8),
        ]
        first = plt.legend(handles=legend1, title=color_col, bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.gca().add_artist(first)
        plt.legend(handles=legend2, title='marker', bbox_to_anchor=(1.02, 0.55), loc='upper left')
    else:
        plt.scatter(df[x_col], df[y_col], c=[color_map[v] for v in df[color_col]], s=20, alpha=0.8, linewidths=0)
        handles = [Line2D([0], [0], marker='o', color='w', label=str(cat),
                          markerfacecolor=color_map[cat], markersize=8) for cat in cats]
        plt.legend(handles=handles, title=color_col, bbox_to_anchor=(1.02, 1), loc='upper left')

    plt.title(title)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def main():
    args = parse_args()
    input_path = Path(args.input_path)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    png_sister_dir = Path(args.png_sister_dir) if args.png_sister_dir else None
    k_values = [int(x) for x in args.k_values.split(",") if x.strip()]

    df, source_files = load_input(input_path)
    if df["id"].duplicated().any():
        dupes = df.loc[df["id"].duplicated(), "id"].tolist()[:10]
        raise ValueError(f"Duplicate file-derived IDs detected. Example duplicates: {dupes}")

    df["topology_len"] = df["topology_string"].str.len()
    df["cropped_topology"] = df["topology_string"].apply(crop_topology)
    df["cropped_len"] = df["cropped_topology"].str.len()
    df["crop_removed_n"] = df["topology_len"] - df["cropped_len"]
    df["is_empty_after_crop"] = df["cropped_len"] == 0

    failed_ids = detect_failed_ids(png_sister_dir)
    df["is_failed_overlay"] = df["id"].isin(failed_ids)

    feature_dicts = [seq_to_features(seq, ks=(2, 3)) for seq in df["cropped_topology"]]
    vec = DictVectorizer(sparse=False)
    X = vec.fit_transform(feature_dicts)
    feature_names = np.array(vec.get_feature_names_out())
    X = normalize(X, norm="l2")

    reducer = umap.UMAP(
        n_neighbors=args.umap_neighbors,
        min_dist=args.umap_min_dist,
        metric="cosine",
        random_state=args.random_state,
    )
    X_umap = reducer.fit_transform(X)
    df["UMAP1"] = X_umap[:, 0]
    df["UMAP2"] = X_umap[:, 1]

    metrics_rows = []
    cluster_summary_frames = []
    top_feature_frames = []

    for k in k_values:
        km = KMeans(n_clusters=k, random_state=args.random_state, n_init=20)
        km_labels = km.fit_predict(X)
        km_col = f"kmeans_k{k}"
        df[km_col] = km_labels
        km_sil = float(silhouette_score(X, km_labels, metric="cosine")) if len(set(km_labels)) > 1 else np.nan
        metrics_rows.append({"method": km_col, "parameter": k, "silhouette_cosine": km_sil})
        cluster_summary_frames.append(summarize_clusters(df, km_col))
        tf = top_features_per_cluster(X, feature_names, km_labels)
        tf.insert(0, "cluster_method", km_col)
        top_feature_frames.append(tf)
        make_umap_plot(df, "UMAP1", "UMAP2", km_col,
                       f"UMAP colored by {km_col}", outdir / f"umap_{km_col}.png", failed_overlay=False)
        make_umap_plot(df, "UMAP1", "UMAP2", km_col,
                       f"UMAP colored by {km_col} with failed overlay", outdir / f"umap_{km_col}_failed_overlay.png", failed_overlay=True)

        ag = AgglomerativeClustering(n_clusters=k, metric="cosine", linkage="average")
        ag_labels = ag.fit_predict(X)
        ag_col = f"agglo_k{k}"
        df[ag_col] = ag_labels
        ag_sil = float(silhouette_score(X, ag_labels, metric="cosine")) if len(set(ag_labels)) > 1 else np.nan
        metrics_rows.append({"method": ag_col, "parameter": k, "silhouette_cosine": ag_sil})
        cluster_summary_frames.append(summarize_clusters(df, ag_col))
        tf = top_features_per_cluster(X, feature_names, ag_labels)
        tf.insert(0, "cluster_method", ag_col)
        top_feature_frames.append(tf)
        make_umap_plot(df, "UMAP1", "UMAP2", ag_col,
                       f"UMAP colored by {ag_col}", outdir / f"umap_{ag_col}.png", failed_overlay=False)
        make_umap_plot(df, "UMAP1", "UMAP2", ag_col,
                       f"UMAP colored by {ag_col} with failed overlay", outdir / f"umap_{ag_col}_failed_overlay.png", failed_overlay=True)

    db = DBSCAN(metric="cosine", eps=0.12, min_samples=10)
    db_labels = db.fit_predict(X)
    df["dbscan"] = db_labels
    unique_db = set(db_labels)
    db_sil = float(silhouette_score(X, db_labels, metric="cosine")) if len(unique_db - {-1}) > 1 else np.nan
    metrics_rows.append({"method": "dbscan", "parameter": "eps=0.12,min_samples=10", "silhouette_cosine": db_sil})
    cluster_summary_frames.append(summarize_clusters(df, "dbscan"))
    tf = top_features_per_cluster(X, feature_names, db_labels)
    tf.insert(0, "cluster_method", "dbscan")
    top_feature_frames.append(tf)
    make_umap_plot(df, "UMAP1", "UMAP2", "dbscan", "UMAP colored by dbscan", outdir / "umap_dbscan.png", failed_overlay=False)
    make_umap_plot(df, "UMAP1", "UMAP2", "dbscan", "UMAP colored by dbscan with failed overlay", outdir / "umap_dbscan_failed_overlay.png", failed_overlay=True)

    plt.figure(figsize=(9, 6))
    plt.hist(df["topology_len"], bins=50, alpha=0.6, label="raw topology length")
    plt.hist(df["cropped_len"], bins=50, alpha=0.6, label="cropped topology length")
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.title("Raw vs cropped topology lengths")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "length_distributions.png", dpi=300)
    plt.close()

    plt.figure(figsize=(9, 6))
    colors = np.where(df["is_failed_overlay"], "#d95f02", "#1b9e77")
    plt.scatter(df["topology_len"], df["cropped_len"], c=colors, alpha=0.75, s=20)
    plt.xlabel("Raw topology length")
    plt.ylabel("Cropped topology length")
    plt.title("Raw vs cropped lengths")
    handles = [
        Line2D([0], [0], marker='o', color='w', label='normal', markerfacecolor='#1b9e77', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='failed overlay', markerfacecolor='#d95f02', markersize=8),
    ]
    plt.legend(handles=handles)
    plt.tight_layout()
    plt.savefig(outdir / "raw_vs_cropped_lengths_scatter.png", dpi=300)
    plt.close()

    n_for_dendro = min(len(df), args.max_dendrogram)
    rng = np.random.default_rng(args.random_state)
    if len(df) > n_for_dendro:
        idx = np.sort(rng.choice(len(df), size=n_for_dendro, replace=False))
        X_d = X[idx]
    else:
        idx = np.arange(len(df))
        X_d = X
    Z = linkage(X_d, method="average", metric="cosine")
    plt.figure(figsize=(14, 6))
    dendrogram(Z, no_labels=True, color_threshold=None)
    plt.title(f"Hierarchical clustering dendrogram (n={len(idx)})")
    plt.xlabel("Sequences")
    plt.ylabel("Cosine distance")
    plt.tight_layout()
    plt.savefig(outdir / "dendrogram_subsample.png", dpi=300)
    plt.close()

    pd.DataFrame(metrics_rows).to_csv(outdir / "clustering_metrics.csv", index=False)
    pd.concat(cluster_summary_frames, ignore_index=True).to_csv(outdir / "cluster_summaries.csv", index=False)
    pd.concat(top_feature_frames, ignore_index=True).to_csv(outdir / "cluster_top_features.csv", index=False)
    df.to_csv(outdir / "sequence_assignments.csv", index=False)

    dataset_summary = {
        "input_path": str(input_path),
        "n_input_files": int(len(source_files)),
        "n_sequences": int(len(df)),
        "n_failed_overlay": int(df["is_failed_overlay"].sum()),
        "fraction_failed_overlay": float(df["is_failed_overlay"].mean()) if len(df) else 0.0,
        "mean_raw_length": float(df["topology_len"].mean()) if len(df) else 0.0,
        "median_raw_length": float(df["topology_len"].median()) if len(df) else 0.0,
        "mean_cropped_length": float(df["cropped_len"].mean()) if len(df) else 0.0,
        "median_cropped_length": float(df["cropped_len"].median()) if len(df) else 0.0,
        "n_empty_after_crop": int(df["is_empty_after_crop"].sum()),
        "feature_count": int(X.shape[1]),
        "umap_neighbors": int(args.umap_neighbors),
        "umap_min_dist": float(args.umap_min_dist),
        "k_values": ",".join(map(str, k_values)),
    }
    pd.DataFrame([dataset_summary]).to_csv(outdir / "dataset_summary.csv", index=False)

    readme = outdir / "README.txt"
    readme.write_text(
        "Unsupervised topology analysis outputs\n"
        "====================================\n\n"
        "Accepted input:\n"
        "- a directory of one-sequence .3line files\n"
        "- or one multi-sequence .3line file\n\n"
        "Main files:\n"
        "- dataset_summary.csv: overall dataset stats after parsing and cropping\n"
        "- sequence_assignments.csv: one row per sequence with crop results, UMAP coordinates, and cluster IDs\n"
        "- clustering_metrics.csv: silhouette scores for clustering runs\n"
        "- cluster_summaries.csv: cluster sizes, mean lengths, and failed-overlay enrichment\n"
        "- cluster_top_features.csv: most enriched k-mer features per cluster\n"
        "- umap_*.png: 2D UMAP visualizations for each clustering run\n"
        "- *_failed_overlay.png: same UMAP plots with failed PNG-linked sequences marked by x\n"
        "- dendrogram_subsample.png: hierarchical clustering dendrogram on all or a reproducible subsample\n"
        "- length_distributions.png, raw_vs_cropped_lengths_scatter.png: crop-effect diagnostics\n\n"
        "Crop rule:\n"
        "- remove leading S\n"
        "- then remove leading O\n"
        "- then remove trailing I until the last character is not I\n"
    )

    print(f"Done. Results written to: {outdir}")


if __name__ == "__main__":
    main()
