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

try:
    import gseapy as gp
except Exception:
    gp = None


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


def slim_feature(feature: str) -> str:
    if pd.isna(feature):
        return "Unknown"

    value = str(feature).strip()
    value_lower = value.lower()

    if value in {"5' UTR", "3' UTR"}:
        return "UTR"
    if value in {"Exon", "Intron"}:
        return "Gene Body"
    if value == "Intergenic" or value == "Distal Intergenic" or value_lower.startswith("downstream"):
        return "Intergenic"

    return value


def sample_id_from_file(path: str) -> str:
    base = os.path.basename(path)
    # output from HOMER_ANNOTATEPEAKS: <sample>.annotatePeaks.txt
    return re.sub(r"\.annotatePeaks\.txt$", "", base)


def condition_id_from_sample(sample: str) -> str:
    """Collapse replicate sample IDs to condition IDs (e.g. X_R1 -> X)."""
    if pd.isna(sample):
        return "Unknown"
    return re.sub(r"_R\d+$", "", str(sample))


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


def clean_gene_tokens(raw_value: str) -> List[str]:
    if pd.isna(raw_value):
        return []

    # HOMER can separate genes with commas/semicolons/slashes and include extra labels.
    text = str(raw_value).strip()
    tokens = re.split(r"[,;/|]+", text)
    cleaned = []
    for token in tokens:
        item = token.strip()
        item = re.sub(r"\s+", " ", item)
        if not item:
            continue

        item_lower = item.lower()
        if item_lower in {"na", "n/a", "none", "unknown", "-", "."}:
            continue

        # Keep common gene symbol characters and avoid trailing metadata noise.
        item = re.sub(r"[^A-Za-z0-9_.\-]", "", item)
        if not item:
            continue

        cleaned.append(item)

    return cleaned


def extract_gene_list(df: pd.DataFrame) -> List[str]:
    gene_col = pick_column(
        list(df.columns),
        [
            "Gene Name",
            "GeneName",
            "gene_name",
            "Nearest PromoterID",
            "Gene ID",
            "gene id",
            "Symbol",
            "gene",
        ],
    )

    if not gene_col:
        return []

    genes: List[str] = []
    for value in df[gene_col].dropna().tolist():
        genes.extend(clean_gene_tokens(value))

    return sorted(set(genes))


def parse_homer_table(path: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str], Dict[str, float]]:
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)

    if df.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), [], {"n_peaks": 0, "mean_gc": float("nan")}

    annotation_col = pick_column(list(df.columns), ["Annotation"])
    if not annotation_col:
        raise ValueError(f"Could not find 'Annotation' column in {path}")

    gc_col = pick_column(list(df.columns), ["GC%", "gc", "gc content"])

    df["feature"] = df[annotation_col].map(canonical_feature)

    n_peaks = int(df.shape[0])
    mean_gc = float("nan")

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

    gc_per_peak = pd.DataFrame(columns=["gc_percent"])
    if gc_col:
        gc_series = pd.to_numeric(df[gc_col], errors="coerce")
        # HOMER often stores GC as percentages already (e.g. 42.1).
        # If values are mostly <= 1, treat them as fractions and convert.
        if gc_series.notna().any() and gc_series.max(skipna=True) <= 1.0:
            gc_series = gc_series * 100.0

        mean_gc = gc_series.mean()

        gc_per_peak = pd.DataFrame({"gc_percent": gc_series}).dropna()

    genes = extract_gene_list(df)

    stats = {
        "n_peaks": n_peaks,
        "mean_gc": mean_gc,
    }
    return raw, feature, gc_per_peak, genes, stats


def run_functional_enrichment(gene_sets: Dict[str, List[str]]) -> pd.DataFrame:
    columns = [
        "sample",
        "term",
        "adjusted_p_value",
        "p_value",
        "odds_ratio",
        "combined_score",
        "overlap",
        "genes",
        "gene_count",
    ]

    if gp is None:
        return pd.DataFrame(columns=columns)

    rows = []
    libraries = ["GO_Biological_Process_2023"]

    for sample, genes in sorted(gene_sets.items()):
        unique_genes = sorted(set(genes))
        if len(unique_genes) < 5:
            continue

        try:
            enr = gp.enrichr(
                gene_list=unique_genes,
                gene_sets=libraries,
                organism="Human",
                outdir=None,
                no_plot=True,
            )
        except Exception:
            continue

        if enr is None or not hasattr(enr, "results") or enr.results is None or enr.results.empty:
            continue

        df = enr.results.copy()
        if "Adjusted P-value" in df.columns:
            df = df.sort_values("Adjusted P-value", ascending=True)
        top_df = df.head(20)

        for _, rec in top_df.iterrows():
            gene_tokens = [g.strip() for g in str(rec.get("Genes", "")).split(";") if g.strip()]
            rows.append(
                {
                    "sample": sample,
                    "term": rec.get("Term", ""),
                    "adjusted_p_value": rec.get("Adjusted P-value", float("nan")),
                    "p_value": rec.get("P-value", float("nan")),
                    "odds_ratio": rec.get("Odds Ratio", float("nan")),
                    "combined_score": rec.get("Combined Score", float("nan")),
                    "overlap": rec.get("Overlap", ""),
                    "genes": rec.get("Genes", ""),
                    "gene_count": len(gene_tokens),
                }
            )

    if not rows:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(rows, columns=columns)


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
    gc_rows = []
    gene_sets: Dict[str, List[str]] = {}
    sample_stats = []

    for input_path in input_files:
        sample_id = sample_id_from_file(input_path)
        raw, feature, gc_per_peak, genes, stats = parse_homer_table(input_path)

        if not raw.empty:
            raw.insert(0, "sample", sample_id)
            raw_rows.append(raw)

        if not feature.empty:
            feature.insert(0, "sample", sample_id)
            feature_rows.append(feature)

        if not gc_per_peak.empty:
            gc_per_peak.insert(0, "sample", sample_id)
            gc_rows.append(gc_per_peak)

        gene_sets[sample_id] = genes

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
    gc_peak_df = pd.concat(gc_rows, ignore_index=True) if gc_rows else pd.DataFrame(columns=["sample", "gc_percent"])
    stats_df = pd.DataFrame(sample_stats).sort_values("sample")

    raw_out = f"{args.out_prefix}.raw_annotation_summary.tsv"
    feature_out = f"{args.out_prefix}.feature_summary.tsv"
    stats_out = f"{args.out_prefix}.sample_stats.tsv"
    gc_peak_out = f"{args.out_prefix}.gc_per_peak.tsv"
    gc_sample_out = f"{args.out_prefix}.gc_by_sample.tsv"
    gc_condition_out = f"{args.out_prefix}.gc_by_condition.tsv"
    enrichment_out = f"{args.out_prefix}.functional_enrichment.tsv"

    raw_df.to_csv(raw_out, sep="\t", index=False)
    feature_df.to_csv(feature_out, sep="\t", index=False)
    stats_df.to_csv(stats_out, sep="\t", index=False)
    gc_peak_df.to_csv(gc_peak_out, sep="\t", index=False)

    if gc_peak_df.empty:
        gc_summary_df = pd.DataFrame(columns=["sample", "n_peaks_with_gc", "mean_gc_percent", "median_gc_percent", "sd_gc_percent"])
    else:
        gc_summary_df = (
            gc_peak_df.groupby("sample", as_index=False)
            .agg(
                n_peaks_with_gc=("gc_percent", "count"),
                mean_gc_percent=("gc_percent", "mean"),
                median_gc_percent=("gc_percent", "median"),
                sd_gc_percent=("gc_percent", "std"),
            )
            .sort_values("sample")
        )
    gc_summary_df.to_csv(gc_sample_out, sep="\t", index=False)

    if gc_peak_df.empty:
        gc_condition_df = pd.DataFrame(
            columns=["condition", "n_peaks_with_gc", "mean_gc_percent", "median_gc_percent", "sd_gc_percent"]
        )
    else:
        gc_condition_df = gc_peak_df.copy()
        gc_condition_df["condition"] = gc_condition_df["sample"].map(condition_id_from_sample)
        gc_condition_df = (
            gc_condition_df.groupby("condition", as_index=False)
            .agg(
                n_peaks_with_gc=("gc_percent", "count"),
                mean_gc_percent=("gc_percent", "mean"),
                median_gc_percent=("gc_percent", "median"),
                sd_gc_percent=("gc_percent", "std"),
            )
            .sort_values("condition")
        )
    gc_condition_df.to_csv(gc_condition_out, sep="\t", index=False)

    enrichment_df = run_functional_enrichment(gene_sets)
    enrichment_df.to_csv(enrichment_out, sep="\t", index=False)

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
    plt.close(fig)

    # Slim version: collapse related features and use half-width figure.
    slim_pivot = pivot.copy()
    slim_pivot.columns = [slim_feature(col) for col in slim_pivot.columns]
    slim_pivot = slim_pivot.T.groupby(level=0).sum().T

    slim_preferred_order = [
        "Promoter",
        "UTR",
        "Gene Body",
        "TTS",
        "Non-coding",
        "Intergenic",
        "Unknown",
    ]
    slim_observed = slim_pivot.columns.tolist()
    slim_ordered_features = [f for f in slim_preferred_order if f in slim_observed]
    slim_ordered_features.extend(sorted([f for f in slim_observed if f not in slim_ordered_features]))
    slim_pivot = slim_pivot.reindex(columns=slim_ordered_features, fill_value=0.0)

    slim_percent_table_out = f"{args.out_prefix}.feature_percent_table.slim_version.tsv"
    slim_pivot.reset_index().to_csv(slim_percent_table_out, sep="\t", index=False)

    slim_default_colors = {
        "Promoter": "#F28E2B",
        "UTR": "#4E79A7",
        "Gene Body": "#59A14F",
        "TTS": "#E15759",
        "Non-coding": "#76B7B2",
        "Intergenic": "#BAB0AC",
        "Unknown": "#9D9D9D",
    }
    slim_colors = [slim_default_colors.get(feature, "#808080") for feature in slim_pivot.columns]

    slim_fig_w = max(4, fig_w / 2)
    fig, ax = plt.subplots(figsize=(slim_fig_w, 6), dpi=150)

    slim_pivot.plot(kind="bar", stacked=True, ax=ax, color=slim_colors, edgecolor="white", linewidth=0.4)

    ax.set_ylabel("Binding (% peaks)")
    ax.set_xlabel("Replicate")
    ax.set_ylim(0, 116)
    ax.set_title("Peak Genomic Feature Composition (Slim Version)")
    ax.legend(title="Feature", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    for i, sample in enumerate(slim_pivot.index):
        n_peaks = int(stats_lookup.loc[sample, "n_peaks"]) if sample in stats_lookup.index else 0
        mean_gc = stats_lookup.loc[sample, "mean_gc_percent"] if sample in stats_lookup.index else float("nan")

        if pd.notna(mean_gc):
            label = f"N={n_peaks}\nGC={mean_gc:.1f}%"
        else:
            label = f"N={n_peaks}"

        ax.text(i, 102.0, label, ha="center", va="bottom", fontsize=8, fontweight="bold")

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    slim_png_out = f"{args.out_prefix}.stacked_bar.slim_version.png"
    slim_pdf_out = f"{args.out_prefix}.stacked_bar.slim_version.pdf"

    plt.savefig(slim_png_out, dpi=200)
    plt.savefig(slim_pdf_out)
    plt.close(fig)

    # GC content visualization by sample.
    if not gc_peak_df.empty:
        gc_plot_png = f"{args.out_prefix}.gc_by_sample.png"
        gc_plot_pdf = f"{args.out_prefix}.gc_by_sample.pdf"
        gc_violin_png = f"{args.out_prefix}.gc_by_sample_violin.png"
        gc_violin_pdf = f"{args.out_prefix}.gc_by_sample_violin.pdf"

        sample_order = sorted(gc_peak_df["sample"].dropna().unique().tolist())
        gc_vectors = [gc_peak_df.loc[gc_peak_df["sample"] == s, "gc_percent"].dropna().tolist() for s in sample_order]

        fig_w = max(8, 0.8 * len(sample_order) + 4)
        fig, ax = plt.subplots(figsize=(fig_w, 6), dpi=150)
        bp = ax.boxplot(gc_vectors, patch_artist=True, showfliers=False)
        for box in bp["boxes"]:
            box.set(facecolor="#4E79A7", alpha=0.65, edgecolor="#2F2F2F")

        ax.set_xticks(range(1, len(sample_order) + 1))
        ax.set_xticklabels(sample_order, rotation=45, ha="right")
        ax.set_ylabel("GC content (%)")
        ax.set_xlabel("Sample")
        ax.set_title("GC Content Distribution by Sample")
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        plt.tight_layout()
        plt.savefig(gc_plot_png, dpi=200)
        plt.savefig(gc_plot_pdf)
        plt.close(fig)

        # Violin version of GC-by-sample distribution.
        fig, ax = plt.subplots(figsize=(fig_w, 6), dpi=150)
        vp = ax.violinplot(gc_vectors, showmeans=False, showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor("#4E79A7")
            body.set_edgecolor("#2F2F2F")
            body.set_alpha(0.5)
        if "cmedians" in vp:
            vp["cmedians"].set_color("#1F1F1F")
            vp["cmedians"].set_linewidth(1.2)

        ax.set_xticks(range(1, len(sample_order) + 1))
        ax.set_xticklabels(sample_order, rotation=45, ha="right")
        ax.set_ylabel("GC content (%)")
        ax.set_xlabel("Sample")
        ax.set_title("GC Content Distribution by Sample (Violin)")
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        plt.tight_layout()
        plt.savefig(gc_violin_png, dpi=200)
        plt.savefig(gc_violin_pdf)
        plt.close(fig)

        # GC content visualization by condition.
        gc_by_condition = gc_peak_df.copy()
        gc_by_condition["condition"] = gc_by_condition["sample"].map(condition_id_from_sample)
        condition_order = sorted(gc_by_condition["condition"].dropna().unique().tolist())
        condition_vectors = [
            gc_by_condition.loc[gc_by_condition["condition"] == c, "gc_percent"].dropna().tolist()
            for c in condition_order
        ]

        cond_fig_w = max(8, 0.8 * len(condition_order) + 4)
        fig, ax = plt.subplots(figsize=(cond_fig_w, 6), dpi=150)
        bp = ax.boxplot(condition_vectors, patch_artist=True, showfliers=False)
        for box in bp["boxes"]:
            box.set(facecolor="#59A14F", alpha=0.65, edgecolor="#2F2F2F")

        ax.set_xticks(range(1, len(condition_order) + 1))
        ax.set_xticklabels(condition_order, rotation=45, ha="right")
        ax.set_ylabel("GC content (%)")
        ax.set_xlabel("Condition")
        ax.set_title("GC Content Distribution by Condition")
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.gc_by_condition.png", dpi=200)
        plt.savefig(f"{args.out_prefix}.gc_by_condition.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(cond_fig_w, 6), dpi=150)
        vp = ax.violinplot(condition_vectors, showmeans=False, showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor("#59A14F")
            body.set_edgecolor("#2F2F2F")
            body.set_alpha(0.5)
        if "cmedians" in vp:
            vp["cmedians"].set_color("#1F1F1F")
            vp["cmedians"].set_linewidth(1.2)

        ax.set_xticks(range(1, len(condition_order) + 1))
        ax.set_xticklabels(condition_order, rotation=45, ha="right")
        ax.set_ylabel("GC content (%)")
        ax.set_xlabel("Condition")
        ax.set_title("GC Content Distribution by Condition (Violin)")
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.gc_by_condition_violin.png", dpi=200)
        plt.savefig(f"{args.out_prefix}.gc_by_condition_violin.pdf")
        plt.close(fig)
    else:
        # Emit placeholder plots so outputs are always present.
        fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
        ax.text(0.5, 0.5, "No GC column detected in annotatePeaks inputs", ha="center", va="center")
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.gc_by_sample.png", dpi=200)
        plt.savefig(f"{args.out_prefix}.gc_by_sample.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
        ax.text(0.5, 0.5, "No GC column detected in annotatePeaks inputs", ha="center", va="center")
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.gc_by_sample_violin.png", dpi=200)
        plt.savefig(f"{args.out_prefix}.gc_by_sample_violin.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
        ax.text(0.5, 0.5, "No GC column detected in annotatePeaks inputs", ha="center", va="center")
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.gc_by_condition.png", dpi=200)
        plt.savefig(f"{args.out_prefix}.gc_by_condition.pdf")
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
        ax.text(0.5, 0.5, "No GC column detected in annotatePeaks inputs", ha="center", va="center")
        ax.set_axis_off()
        plt.tight_layout()
        plt.savefig(f"{args.out_prefix}.gc_by_condition_violin.png", dpi=200)
        plt.savefig(f"{args.out_prefix}.gc_by_condition_violin.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
