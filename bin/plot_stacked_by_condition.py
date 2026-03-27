#!/usr/bin/env python3

"""Create condition-level stacked bar plots from a wide sample table.

Expected input format:
- One row per sample
- First column is sample name (default: "Sample")
- Remaining numeric columns are categories (e.g., Human, Ad5, No hits)

The script infers condition labels from sample names by default:
  KM_3_2_screen      -> KM
  IgG_control_2_screen -> IgG_control

It outputs:
- <prefix>.condition_means.tsv
- <prefix>.condition_means_collapsed.tsv
- <prefix>.stacked_by_condition.png
- <prefix>.stacked_by_condition.pdf
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a stacked bar plot by condition from wide sample percentages"
    )
    parser.add_argument("--input", required=True, help="Input TSV/CSV file")
    parser.add_argument(
        "--out-prefix",
        default="screen_taxa",
        help="Output file prefix (default: screen_taxa)",
    )
    parser.add_argument(
        "--sample-col",
        default="Sample",
        help="Column containing sample IDs (default: Sample)",
    )
    parser.add_argument(
        "--condition-regex",
        default=r"^(.*?)(?:_\d+)?(?:_\d+)?_screen$",
        help=(
            "Regex with one capture group used to infer condition from sample name. "
            "Default handles names like KM_3_2_screen and IgG_control_2_screen."
        ),
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        help="Optional delimiter override. If omitted, auto-detects tab/comma.",
    )
    parser.add_argument(
        "--other-threshold",
        type=float,
        default=1.0,
        help=(
            "Collapse values below this percent into an 'Other' bucket per condition "
            "(default: 1.0)."
        ),
    )
    return parser.parse_args()


def infer_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline()
    return "\t" if header.count("\t") >= header.count(",") else ","


def infer_condition(sample: str, pattern: re.Pattern[str]) -> str:
    sample = str(sample)
    m = pattern.match(sample)
    if m:
        return m.group(1)
    return sample


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    delim = args.delimiter if args.delimiter is not None else infer_delimiter(input_path)
    df = pd.read_csv(input_path, sep=delim)

    if args.sample_col not in df.columns:
        raise ValueError(
            f"Sample column '{args.sample_col}' not found. Available columns: {list(df.columns)}"
        )

    numeric_cols = [c for c in df.columns if c != args.sample_col]
    if not numeric_cols:
        raise ValueError("No category columns found.")

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    pattern = re.compile(args.condition_regex)
    df["condition"] = df[args.sample_col].astype(str).map(lambda s: infer_condition(s, pattern))

    condition_means = (
        df.groupby("condition", dropna=False)[numeric_cols]
        .mean()
        .fillna(0.0)
        .sort_index()
    )

    # Collapse tiny components to an 'Other' bucket per condition.
    plot_df = condition_means.copy()
    threshold = max(0.0, float(args.other_threshold))
    low_mask = plot_df < threshold
    other_values = plot_df.where(low_mask, 0.0).sum(axis=1)
    plot_df = plot_df.mask(low_mask, 0.0)
    plot_df = plot_df.loc[:, (plot_df > 0).any(axis=0)]
    if (other_values > 0).any():
        plot_df["Other"] = other_values

    # Keep categories in descending overall contribution so legend/stacking are stable.
    ordered_categories = plot_df.mean(axis=0).sort_values(ascending=False).index.tolist()
    if "Other" in ordered_categories:
        ordered_categories = [c for c in ordered_categories if c != "Other"] + ["Other"]
    plot_df = plot_df[ordered_categories]

    out_prefix = Path(args.out_prefix)
    condition_table = out_prefix.with_suffix("")
    condition_means.reset_index().to_csv(
        f"{condition_table}.condition_means.tsv", sep="\t", index=False
    )
    plot_df.reset_index().to_csv(
        f"{condition_table}.condition_means_collapsed.tsv", sep="\t", index=False
    )

    fig_w = max(8, 0.8 * len(plot_df.index) + 4)
    fig, ax = plt.subplots(figsize=(fig_w, 6), dpi=150)

    cmap = plt.cm.get_cmap("tab20", len(ordered_categories))
    colors = [cmap(i) for i in range(len(ordered_categories))]

    plot_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        edgecolor="white",
        linewidth=0.35,
    )

    ax.set_title(f"Stacked Composition by Condition (Other < {threshold:.2f}%)")
    ax.set_xlabel("Condition")
    ax.set_ylabel("Percent")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", linestyle=":", alpha=0.3)
    ax.legend(title="Category", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(f"{condition_table}.stacked_by_condition.png", dpi=200)
    plt.savefig(f"{condition_table}.stacked_by_condition.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
