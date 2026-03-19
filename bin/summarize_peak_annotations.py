#!/usr/bin/env python3

# MIT License
#
# Copyright (c) 2026 Katharina Hayer
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
import glob
import os
import re
from typing import Dict, List, Tuple

import pandas as pd
import matplotlib.pyplot as plt


def canonical_feature(annotation: str) -> str:
    if pd.isna(annotation):
        return "Unknown"

    ann = str(annotation).strip()
    ann_lower = ann.lower()

    if "5' utr" in ann_lower or "5 utr" in ann_lower:
        return "5' UTR"
    if "3' utr" in ann_lower or "3 utr" in ann_lower:
        return "3' UTR"
    if ann_lower.startswith("promoter") or "tss" in ann_lower:
        return "Promoter"
    if ann_lower.startswith("exon"):
        return "Exon"
    if ann_lower.startswith("intron"):
        return "Intron"
    if "intergenic" in ann_lower:
        return "Intergenic"
    if ann_lower.startswith("tts"):
        return "TTS"
    if "non-coding" in ann_lower:
        return "Non-coding"

    # Keep HOMER wording if it does not match a known bucket.
    return ann


def sample_id_from_file(path: str) -> str:
    base = os.path.basename(path)
    # output from HOMER_ANNOTATEPEAKS: <sample>.annotatePeaks.txt
    return re.sub(r"\.annotatePeaks\.txt$", "", base)


def pick_column(columns: List[str], candidates: List[str]) -> str:
    lower_to_original = {c.lower(): c for c in columns}
    for candidate in candidates:
        if candidate.lower() in lower_to_original:
            return lower_to_original[candidate.lower()]
    for c in columns:
        c_lower = c.lower()
        if any(x in c_lower for x in candidates):
            return c
    return ""


def parse_homer_table(path: str) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)

    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), {"n_peaks": 0, "mean_gc": float("nan")}

    annotation_col = pick_column(list(df.columns), ["Annotation"])
    if not annotation_col:
        raise ValueError(f"Could not find 'Annotation' column in {path}")

    gc_col = pick_column(list(df.columns), ["GC%", "gc", "gc content"])

    df["feature"] = df[annotation_col].map(canonical_feature)

    n_peaks = int(df.shape[0])
    mean_gc = pd.to_numeric(df[gc_col], errors="coerce").mean() if gc_col else float("nan")

    # Preserve raw HOMER annotation classes for complete reporting.
    raw = (
        df[annotation_col]
        .fillna("Unknown")
        .value_counts(dropna=False)
        .rename_axis("raw_annotation")
        .reset_index(name="count")
    )
    raw["percent"] = (raw["count"] / n_peaks * 100.0) if n_peaks else 0.0

    # Feature buckets used for stacked bar plot.
    feature = (
        df["feature"]
        .fillna("Unknown")
        .value_counts(dropna=False)
        .rename_axis("feature")
        .reset_index(name="count")
    )
    feature["percent"] = (feature["count"] / n_peaks * 100.0) if n_peaks else 0.0

    stats = {
        "n_peaks": n_peaks,
        "mean_gc": mean_gc,
    }
    return raw, feature, stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize and plot HOMER annotatePeaks results per replicate."
    )
    parser.add_argument(
        "--inputs",
        required=True,
        nargs="+",
        help="One or more HOMER annotatePeaks output files",
    )
    parser.add_argument(
        "--out-prefix",
        default="homer_peak_annotation",
        help="Output file prefix",
    )
    args = parser.parse_args()

    input_files: List[str] = []
    for pattern in args.inputs:
        expanded = glob.glob(pattern)
        if expanded:
            input_files.extend(expanded)
        elif os.path.exists(pattern):
            input_files.append(pattern)

    input_files = sorted(set(input_files))
    if not input_files:
        raise SystemExit("No annotatePeaks input files were found.")

    raw_rows = []
    feature_rows = []
    sample_stats = []

    for input_path in input_files:
        sample_id = sample_id_from_file(input_path)
        raw, feature, stats = parse_homer_table(input_path)

        if not raw.empty:
            raw.insert(0, "sample", sample_id)
            raw_rows.append(raw)

        if not feature.empty:
            feature.insert(0, "sample", sample_id)
            feature_rows.append(feature)

        sample_stats.append(
            {
                "sample": sample_id,
                "n_peaks": stats["n_peaks"],
                "mean_gc_percent": stats["mean_gc"],
                "source_file": os.path.basename(input_path),
            }
        )

    raw_df = pd.concat(raw_rows, ignore_index=True) if raw_rows else pd.DataFrame()
    feature_df = pd.concat(feature_rows, ignore_index=True) if feature_rows else pd.DataFrame()
    stats_df = pd.DataFrame(sample_stats).sort_values("sample")

    raw_out = f"{args.out_prefix}.raw_annotation_summary.tsv"
    feature_out = f"{args.out_prefix}.feature_summary.tsv"
    stats_out = f"{args.out_prefix}.sample_stats.tsv"

    raw_df.to_csv(raw_out, sep="\t", index=False)
    feature_df.to_csv(feature_out, sep="\t", index=False)
    stats_df.to_csv(stats_out, sep="\t", index=False)

    if feature_df.empty:
        raise SystemExit("No feature annotations were parsed from annotatePeaks outputs.")

    # Feature order with explicit UTR split first, then other frequent classes.
    preferred_order = [
        "Promoter",
        "5' UTR",
        "3' UTR",
        "Exon",
        "Intron",
        "TTS",
        "Non-coding",
        "Intergenic",
        "Unknown",
    ]

    observed = feature_df["feature"].dropna().unique().tolist()
    ordered_features = [f for f in preferred_order if f in observed]
    ordered_features.extend(sorted([f for f in observed if f not in ordered_features]))

    pivot = (
        feature_df.pivot_table(index="sample", columns="feature", values="percent", aggfunc="sum", fill_value=0.0)
        .reindex(columns=ordered_features, fill_value=0.0)
        .sort_index()
    )

    # Keep a numeric table used by downstream users.
    percent_table_out = f"{args.out_prefix}.feature_percent_table.tsv"
    pivot.reset_index().to_csv(percent_table_out, sep="\t", index=False)

    # Color map (stable, readable palette).
    default_colors = {
        "Promoter": "#F28E2B",
        "5' UTR": "#86BCB6",
        "3' UTR": "#4E79A7",
        "Exon": "#59A14F",
        "Intron": "#B07AA1",
        "TTS": "#E15759",
        "Non-coding": "#76B7B2",
        "Intergenic": "#BAB0AC",
        "Unknown": "#9D9D9D",
    }
    colors = [default_colors.get(feature, "#808080") for feature in pivot.columns]

    fig_w = max(8, 0.8 * len(pivot.index) + 4)
    fig, ax = plt.subplots(figsize=(fig_w, 6), dpi=150)

    pivot.plot(kind="bar", stacked=True, ax=ax, color=colors, edgecolor="white", linewidth=0.4)

    ax.set_ylabel("Binding (% peaks)")
    ax.set_xlabel("Replicate")
    ax.set_ylim(0, 116)
    ax.set_title("Peak Genomic Feature Composition (HOMER annotatePeaks)")
    ax.legend(title="Feature", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    stats_lookup = stats_df.set_index("sample")
    for i, sample in enumerate(pivot.index):
        n_peaks = int(stats_lookup.loc[sample, "n_peaks"]) if sample in stats_lookup.index else 0
        mean_gc = stats_lookup.loc[sample, "mean_gc_percent"] if sample in stats_lookup.index else float("nan")

        if pd.notna(mean_gc):
            label = f"N={n_peaks}\nGC={mean_gc:.1f}%"
        else:
            label = f"N={n_peaks}"

        ax.text(i, 102.0, label, ha="center", va="bottom", fontsize=8, fontweight="bold")

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    png_out = f"{args.out_prefix}.stacked_bar.png"
    pdf_out = f"{args.out_prefix}.stacked_bar.pdf"

    plt.savefig(png_out, dpi=200)
    plt.savefig(pdf_out)


if __name__ == "__main__":
    main()
