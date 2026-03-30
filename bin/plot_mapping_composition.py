#!/usr/bin/env python3
"""
Plot mapping composition from MultiQC samtools idxstats data.
Generates stacked bar plots showing percentage of reads mapped to:
- Human genome references
- Ad5 (viral contaminant)
- Other (known non-target organisms)
"""

import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt


def parse_count(x, count_re):
    """Extract mapped read count from [count, length] format."""
    if pd.isna(x):
        return 0
    m = count_re.search(str(x))
    return int(m.group(1)) if m else 0


def main():
    infile = 'multiqc_samtools_idxstats.txt'
    outdir = 'results_human/mapping_summary'
    
    if not os.path.exists(infile):
        print(f"Error: {infile} not found", file=sys.stderr)
        sys.exit(1)
    
    os.makedirs(outdir, exist_ok=True)
    
    # Read data
    df = pd.read_csv(infile, sep='\t')
    count_re = re.compile(r"\[\s*(\d+)\s*,")
    
    # Identify reference columns
    ref_cols = [c for c in df.columns if c != 'Sample']
    non_human_known = {
        'Ad5',
        'HSV_HSJN555585.1',
        'hRSV_A2',
        'Myco_NC_010163.1',
        'Myco_NC_014921.1',
        'Myco_NC_013511.1',
        'Myco_NC_017519.1',
        '*'
    }
    
    ad5_cols = [c for c in ref_cols if c == 'Ad5']
    other_cols = [c for c in ref_cols if c in non_human_known and c not in {'Ad5', '*'}]
    human_cols = [c for c in ref_cols if c not in non_human_known]
    
    # Extract and sum mapped counts by category
    counts = pd.DataFrame({'Sample': df['Sample']})
    counts['Human'] = df[human_cols].applymap(lambda x: parse_count(x, count_re)).sum(axis=1)
    counts['Ad5'] = df[ad5_cols].applymap(lambda x: parse_count(x, count_re)).sum(axis=1) if ad5_cols else 0
    counts['Other'] = df[other_cols].applymap(lambda x: parse_count(x, count_re)).sum(axis=1)
    counts['Mapped_total'] = counts[['Human', 'Ad5', 'Other']].sum(axis=1)
    
    # Calculate percentages
    pct = counts.copy()
    for col in ['Human', 'Ad5', 'Other']:
        pct[col] = (pct[col] / pct['Mapped_total'].replace(0, pd.NA) * 100).fillna(0)
    
    # Save sample-level percentages
    pct[['Sample', 'Human', 'Ad5', 'Other', 'Mapped_total']].to_csv(
        os.path.join(outdir, 'mapped_category_percent_by_sample.tsv'), sep='\t', index=False)
    
    # Define colors
    colors = {'Human': '#4C78A8', 'Ad5': '#F58518', 'Other': '#54A24B'}
    
    # Plot 1: Sample-level stacked bar chart
    fig, ax = plt.subplots(figsize=(11, 6), dpi=150)
    x = range(len(pct))
    bottom = [0] * len(pct)
    for col in ['Human', 'Ad5', 'Other']:
        ax.bar(x, pct[col], bottom=bottom, label=col, color=colors[col], width=0.85)
        bottom = [b + v for b, v in zip(bottom, pct[col])]
    
    ax.set_ylim(0, 100)
    ax.set_ylabel('Mapped reads (%)')
    ax.set_title('Mapped Read Composition by Sample (Ad5 vs Human vs Other)')
    ax.set_xticks(list(x))
    ax.set_xticklabels(pct['Sample'], rotation=45, ha='right')
    ax.legend(frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.15))
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'mapped_category_stacked_by_sample.png'))
    fig.savefig(os.path.join(outdir, 'mapped_category_stacked_by_sample.pdf'))
    plt.close(fig)
    
    # Plot 2: Condition-mean stacked bar chart
    cond = pct.copy()
    cond['Condition'] = cond['Sample'].str.replace(r'_(?:R)?\d+$', '', regex=True)
    cond_mean = cond.groupby('Condition', as_index=False)[['Human', 'Ad5', 'Other']].mean()
    cond_mean.to_csv(os.path.join(outdir, 'mapped_category_percent_by_condition_mean.tsv'), sep='\t', index=False)
    
    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=150)
    x = range(len(cond_mean))
    bottom = [0] * len(cond_mean)
    for col in ['Human', 'Ad5', 'Other']:
        ax.bar(x, cond_mean[col], bottom=bottom, label=col, color=colors[col], width=0.85)
        bottom = [b + v for b, v in zip(bottom, cond_mean[col])]
    
    ax.set_ylim(0, 100)
    ax.set_ylabel('Mapped reads (%)')
    ax.set_title('Mean Mapped Read Composition by Condition (Ad5 vs Human vs Other)')
    ax.set_xticks(list(x))
    ax.set_xticklabels(cond_mean['Condition'], rotation=45, ha='right')
    ax.legend(frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.15))
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'mapped_category_stacked_by_condition_mean.png'))
    fig.savefig(os.path.join(outdir, 'mapped_category_stacked_by_condition_mean.pdf'))
    plt.close(fig)
    
    print(f"Plots saved to {outdir}")
    print(f"Human refs: {len(human_cols)}")
    print(f"Other refs: {other_cols}")


if __name__ == '__main__':
    main()
