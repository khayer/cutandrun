#!/usr/bin/env python3
"""
Create motif comparison tables across conditions.

Author: Katharina Hayer
Co-created with: GitHub Copilot (Claude Sonnet 4.5)

This script parses Homer motif results and creates comparison tables
showing which motifs are found in which conditions with their % Target values.
"""

import os
import sys
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

try:
    import cairosvg
except Exception:
    cairosvg = None


def parse_known_motifs(filepath):
    """Parse Homer knownResults.txt file."""
    motifs = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
        # Skip header
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            
            motif_name = parts[0]
            # Extract just the base motif name (before the first parenthesis)
            base_name = motif_name.split('(')[0]
            consensus = parts[1]
            pvalue = parts[2]
            percent_target = parts[6]  # "X.XX%"
            
            motifs.append({
                'motif_name': base_name,
                'full_name': motif_name,
                'consensus': consensus,
                'pvalue': pvalue,
                'percent_target': float(percent_target.rstrip('%'))
            })
    return motifs


def parse_denovo_motifs(filepath):
    """Parse Homer homerMotifs.all.motifs file."""
    motifs = []
    motif_dir = Path(filepath).parent / 'homerResults'
    
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Parse header line
                # Format: >CONSENSUS	ID	score	...	T:X(Y%),B:...	...
                parts = line.strip().split('\t')
                consensus = parts[0][1:]  # Remove '>'
                motif_id = parts[1] if len(parts) > 1 else consensus
                
                # Find the part with T:X(Y%)
                target_info = None
                for part in parts:
                    if part.startswith('T:'):
                        target_info = part
                        break
                
                if target_info:
                    # Extract percentage from T:X(Y%),B:...
                    match = re.search(r'T:[\d.]+\(([\d.]+)%\)', target_info)
                    if match:
                        percent = float(match.group(1))
                        
                        # Extract p-value if present
                        pvalue = None
                        pvalue_match = re.search(r'P:([\de-]+)', target_info)
                        if pvalue_match:
                            pvalue = pvalue_match.group(1)
                        
                        # Determine SVG path
                        # motif_id is like "1-CONSENSUS" -> extract the number
                        motif_num = motif_id.split('-')[0]
                        svg_path = motif_dir / f'motif{motif_num}.logo.svg'
                        
                        motifs.append({
                            'motif_id': motif_id,
                            'consensus': consensus,
                            'pvalue': pvalue,
                            'percent_target': percent,
                            'svg_path': str(svg_path) if svg_path.exists() else None
                        })
    return motifs


def create_known_motif_table(motif_dir, output_file):
    """Create comparison table for known motifs."""
    
    # Find all condition directories
    conditions = []
    condition_motifs = {}
    
    consensus_dir = Path(motif_dir) / 'consensus_peaks'
    merged_dir = Path(motif_dir) / 'merged_peaks'
    
    # Process consensus peaks
    if consensus_dir.exists():
        for cond_dir in sorted(consensus_dir.iterdir()):
            if cond_dir.is_dir() and cond_dir.name.endswith('_motifs'):
                condition = cond_dir.name.replace('_motifs', '')
                known_file = cond_dir / 'knownResults.txt'
                if known_file.exists():
                    conditions.append(condition)
                    condition_motifs[condition] = parse_known_motifs(known_file)
    
    # Process merged peaks
    if merged_dir.exists():
        merged_motifs_dir = merged_dir / 'merged_peaks_motifs'
        if merged_motifs_dir.exists():
            known_file = merged_motifs_dir / 'knownResults.txt'
            if known_file.exists():
                conditions.append('merged_peaks')
                condition_motifs['merged_peaks'] = parse_known_motifs(known_file)
    
    if not conditions:
        print("No known motif results found")
        return
    
    # Separate consensus conditions from merged_peaks
    consensus_conditions = [c for c in conditions if c != 'merged_peaks']
    
    # Collect all unique motifs by base name and consensus
    all_motifs = {}  # key: (motif_name, consensus), value: full_name example
    for cond, motifs in condition_motifs.items():
        for m in motifs:
            key = (m['motif_name'], m['consensus'])
            if key not in all_motifs:
                all_motifs[key] = m['full_name']
    
    # Create table
    rows = []
    for (motif_name, consensus), full_name in sorted(all_motifs.items()):
        row = {
            'Motif': motif_name,
            'Consensus': consensus,
        }
        
        # Store raw values for calculations
        percent_values = []
        
        # Add values for each condition (separate % and p-value columns)
        for cond in conditions:
            motifs_in_cond = condition_motifs[cond]
            # Find matching motif
            match = None
            for m in motifs_in_cond:
                if m['motif_name'] == motif_name and m['consensus'] == consensus:
                    match = m
                    break
            
            if match:
                row[cond] = f"{match['percent_target']:.2f}%"
                # Mark non-significant p-values as -1
                if match['pvalue'] in ['1e0', '1.0', '1']:
                    row[f"{cond}_pval"] = '-1'
                else:
                    row[f"{cond}_pval"] = match['pvalue']
                # Only include consensus conditions in average calculation
                if cond != 'merged_peaks':
                    percent_values.append(match['percent_target'])
            else:
                row[cond] = '-'
                row[f"{cond}_pval"] = '-'
        
        rows.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(rows)
    
    # Add average % target (excluding merged_peaks and zero values)
    def calc_average(row):
        values = []
        for cond in consensus_conditions:
            if row[cond] != '-':
                # Extract percentage from "XX.XX%" format
                pct_str = row[cond].rstrip('%')
                try:
                    pct = float(pct_str)
                    if pct > 0:
                        values.append(pct)
                except:
                    pass
        return sum(values) / len(values) if values else 0.0
    
    df['Avg_Target'] = df.apply(calc_average, axis=1)
    
    # Add count of conditions with this motif (excluding merged_peaks, "-", and "0.00%")
    def count_conditions(row):
        count = 0
        for cond in consensus_conditions:
            if row[cond] != '-':
                # Extract percentage
                pct_str = row[cond].rstrip('%')
                try:
                    if float(pct_str) > 0:
                        count += 1
                except:
                    pass
        return count
    
    df['Conditions_Found'] = df.apply(count_conditions, axis=1)
    
    # Format average
    df['Avg_Target'] = df['Avg_Target'].apply(lambda x: f"{x:.2f}%")
    
    # Filter out rows where all p-values are -1
    pval_cols = [f"{cond}_pval" for cond in conditions]
    def has_significant_pval(row):
        for col in pval_cols:
            if row[col] not in ['-1', '-']:
                return True
        return False
    
    df = df[df.apply(has_significant_pval, axis=1)]
    
    # Sort by Conditions_Found (descending) then by Avg_Target (descending)
    df['_sort_avg'] = df['Avg_Target'].apply(lambda x: float(x.rstrip('%')))
    df = df.sort_values(['Conditions_Found', '_sort_avg'], ascending=[False, False])
    df = df.drop(columns=['_sort_avg'])
    
    # Reorder columns: Motif, Consensus, then alternating % and pval for each condition, then metrics
    base_cols = ['Motif', 'Consensus']
    condition_cols = []
    for cond in conditions:
        condition_cols.append(cond)
        condition_cols.append(f"{cond}_pval")
    metric_cols = ['Conditions_Found', 'Avg_Target']
    
    cols = base_cols + condition_cols + metric_cols
    df = df[cols]
    
    # Save to file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Created known motif comparison table: {output_file}")
    print(f"  {len(rows)} unique motifs across {len(conditions)} conditions")
    
    return df


def create_denovo_motif_table(motif_dir, output_file):
    """Create comparison table for de novo motifs."""
    
    # Find all condition directories
    conditions = []
    condition_motifs = {}
    
    consensus_dir = Path(motif_dir) / 'consensus_peaks'
    merged_dir = Path(motif_dir) / 'merged_peaks'
    
    # Process consensus peaks
    if consensus_dir.exists():
        for cond_dir in sorted(consensus_dir.iterdir()):
            if cond_dir.is_dir() and cond_dir.name.endswith('_motifs'):
                condition = cond_dir.name.replace('_motifs', '')
                denovo_file = cond_dir / 'homerMotifs.all.motifs'
                if denovo_file.exists():
                    conditions.append(condition)
                    condition_motifs[condition] = parse_denovo_motifs(denovo_file)
    
    # Process merged peaks
    if merged_dir.exists():
        merged_motifs_dir = merged_dir / 'merged_peaks_motifs'
        if merged_motifs_dir.exists():
            denovo_file = merged_motifs_dir / 'homerMotifs.all.motifs'
            if denovo_file.exists():
                conditions.append('merged_peaks')
                condition_motifs['merged_peaks'] = parse_denovo_motifs(denovo_file)
    
    if not conditions:
        print("No de novo motif results found")
        return
    
    # Separate consensus conditions from merged_peaks
    consensus_conditions = [c for c in conditions if c != 'merged_peaks']
    
    # Collect all unique motifs by consensus
    all_consensus = set()
    consensus_to_svg = {}  # Map consensus to SVG paths
    
    for cond, motifs in condition_motifs.items():
        for m in motifs:
            all_consensus.add(m['consensus'])
            # Store the first SVG path found for this consensus
            if m['consensus'] not in consensus_to_svg and m.get('svg_path'):
                consensus_to_svg[m['consensus']] = m['svg_path']
    
    # Create table
    rows = []
    for consensus in sorted(all_consensus):
        row = {
            'Consensus': consensus,
        }
        
        # Add SVG path if available
        row['SVG_Path'] = consensus_to_svg.get(consensus, '-')
        
        # Add values for each condition (separate % and p-value columns)
        for cond in conditions:
            motifs_in_cond = condition_motifs[cond]
            # Find matching motif
            match = None
            for m in motifs_in_cond:
                if m['consensus'] == consensus:
                    match = m
                    break
            
            if match:
                row[cond] = f"{match['percent_target']:.2f}%"
                if match['pvalue']:
                    # Mark non-significant p-values as -1
                    if match['pvalue'] in ['1e0', '1.0', '1']:
                        row[f"{cond}_pval"] = '-1'
                    else:
                        row[f"{cond}_pval"] = match['pvalue']
                else:
                    row[f"{cond}_pval"] = '-'
            else:
                row[cond] = '-'
                row[f"{cond}_pval"] = '-'
        
        rows.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(rows)
    
    # Add average % target (excluding merged_peaks and zero values)
    def calc_average(row):
        values = []
        for cond in consensus_conditions:
            if row[cond] != '-':
                # Extract percentage from "XX.XX%" format
                pct_str = row[cond].rstrip('%')
                try:
                    pct = float(pct_str)
                    if pct > 0:
                        values.append(pct)
                except:
                    pass
        return sum(values) / len(values) if values else 0.0
    
    df['Avg_Target'] = df.apply(calc_average, axis=1)
    
    # Add count of conditions with this motif (excluding merged_peaks, "-", and "0.00%")
    def count_conditions(row):
        count = 0
        for cond in consensus_conditions:
            if row[cond] != '-':
                # Extract percentage
                pct_str = row[cond].rstrip('%')
                try:
                    if float(pct_str) > 0:
                        count += 1
                except:
                    pass
        return count
    
    df['Conditions_Found'] = df.apply(count_conditions, axis=1)
    
    # Format average
    df['Avg_Target'] = df['Avg_Target'].apply(lambda x: f"{x:.2f}%")
    
    # Sort by Conditions_Found (descending) then by Avg_Target (descending)
    df['_sort_avg'] = df['Avg_Target'].apply(lambda x: float(x.rstrip('%')))
    df = df.sort_values(['Conditions_Found', '_sort_avg'], ascending=[False, False])
    df = df.drop(columns=['_sort_avg'])
    
    # Reorder columns: Consensus, then alternating % and pval for each condition, then metrics, then SVG_Path last
    base_cols = ['Consensus']
    condition_cols = []
    for cond in conditions:
        condition_cols.append(cond)
        condition_cols.append(f"{cond}_pval")
    metric_cols = ['Conditions_Found', 'Avg_Target']
    
    cols = base_cols + condition_cols + metric_cols + ['SVG_Path']
    df = df[cols]
    
    # Save to file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Created de novo motif comparison table: {output_file}")
    print(f"  {len(rows)} unique motifs across {len(conditions)} conditions")
    
    return df


def _safe_int(value, default=0):
    try:
        return int(value)
    except Exception:
        return default


def _svg_to_png(svg_file, out_png):
    """Convert SVG to PNG using cairosvg or rsvg-convert fallback."""
    if not svg_file or not Path(svg_file).exists():
        return False

    if cairosvg is not None:
        try:
            cairosvg.svg2png(url=str(svg_file), write_to=str(out_png), output_height=700)
            return Path(out_png).exists()
        except Exception:
            pass

    rsvg_convert = shutil.which('rsvg-convert')
    if not rsvg_convert and Path('/usr/bin/rsvg-convert').exists():
        rsvg_convert = '/usr/bin/rsvg-convert'

    if not rsvg_convert:
        return False
    try:
        subprocess.run(
            [rsvg_convert, '-h', '700', '-o', str(out_png), str(svg_file)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return Path(out_png).exists()
    except Exception:
        return False


def _render_logo_table_pdf(rows, title, output_pdf):
    """Render a table-like PDF with a sequence-logo column using matplotlib axes."""
    n_rows = len(rows)
    n_cols = 7
    col_labels = ['Rank', 'Motif', 'Consensus', 'P-value', '% Targets', '% Background', 'SVG Logo']
    max_motif_len = max([max(1, int(r.get('motif_len', 1))) for r in rows]) if rows else 1

    fig_w = 16.5
    fig_h = max(3.0, 1.15 * (n_rows + 1) + 1.2)
    fig, axes = plt.subplots(
        nrows=n_rows + 1,
        ncols=n_cols,
        figsize=(fig_w, fig_h),
        gridspec_kw={'width_ratios': [0.6, 3.1, 1.5, 1.2, 1.2, 1.4, 4.0]},
    )

    if n_rows == 0:
        axes = [axes]

    # Header
    for c in range(n_cols):
        ax = axes[0][c]
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor('#d9d9d9')
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
            spine.set_color('black')
        ax.text(0.5, 0.5, col_labels[c], ha='center', va='center', fontsize=8, fontweight='bold')

    # Body rows
    for r, row in enumerate(rows, start=1):
        text_vals = [
            str(row['rank']),
            row['motif'],
            row['consensus'],
            row['pval'],
            row['pct_targets'],
            row['pct_bg'],
        ]

        for c in range(n_cols):
            ax = axes[r][c]
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_facecolor('#f7f7f7' if r % 2 == 1 else 'white')
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(0.6)
                spine.set_color('black')

            if c < 6:
                ax.text(0.5, 0.5, text_vals[c], ha='center', va='center', fontsize=7, wrap=True)
            else:
                logo_png = row.get('logo_png')
                if logo_png and Path(logo_png).exists():
                    try:
                        img = plt.imread(logo_png)
                        motif_len = max(1, int(row.get('motif_len', 1)))
                        width_frac = min(1.0, motif_len / max_motif_len)
                        ax.set_xlim(0, 1)
                        ax.set_ylim(0, 1)
                        # Keep logos left-aligned and width-scaled by motif length
                        # so base positions align across rows.
                        ax.imshow(
                            img,
                            interpolation='lanczos',
                            extent=(0, width_frac, 0, 1),
                            aspect='auto',
                        )
                    except Exception:
                        ax.text(0.5, 0.5, 'SVG render failed', ha='center', va='center', fontsize=7)
                else:
                    ax.text(0.5, 0.5, 'No logo', ha='center', va='center', fontsize=7)

    fig.suptitle(title, fontsize=12, fontweight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(output_pdf, format='pdf', dpi=300)
    plt.close(fig)


def create_known_motif_pdf(homer_dir, output_pdf, top_n=5):
    """Create a per-condition PDF table from Homer knownResults.txt output."""
    homer_dir = Path(homer_dir)
    known_file = homer_dir / 'knownResults.txt'

    if not known_file.exists():
        print(f"Skipping PDF (missing knownResults.txt): {homer_dir}")
        return False

    motifs = []
    with open(known_file, 'r') as fh:
        lines = fh.readlines()
        for i, line in enumerate(lines[1:top_n + 1], start=1):
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            svg_file = homer_dir / 'knownResults' / f"known{i}.logo.svg"
            motifs.append({
                'rank': len(motifs) + 1,
                'motif': parts[0],
                'consensus': parts[1],
                'pval': parts[2],
                'pct_targets': parts[6],
                'pct_bg': parts[8],
                'svg_file': str(svg_file) if svg_file.exists() else None,
                'motif_len': len(parts[1]),
            })

    if not motifs:
        print(f"Skipping PDF (no parsed motifs): {known_file}")
        return False

    with tempfile.TemporaryDirectory() as tmp_dir:
        for idx, m in enumerate(motifs, start=1):
            out_png = Path(tmp_dir) / f"known_{idx}.png"
            if _svg_to_png(m.get('svg_file'), out_png):
                m['logo_png'] = str(out_png)
            else:
                m['logo_png'] = None

        _render_logo_table_pdf(
            motifs,
            f"HOMER Known Motifs: {homer_dir.name} (Top {len(motifs)})",
            output_pdf,
        )

    print(f"Created motif table PDF: {output_pdf}")
    return True


def create_denovo_motif_pdf(homer_dir, output_pdf, top_n=5):
    """Create a per-condition PDF table from Homer de novo motif outputs."""
    homer_dir = Path(homer_dir)
    denovo_file = homer_dir / 'homerMotifs.all.motifs'
    homer_results_dir = homer_dir / 'homerResults'

    if not denovo_file.exists():
        print(f"Skipping de novo PDF (missing homerMotifs.all.motifs): {homer_dir}")
        return False

    motifs = []
    with open(denovo_file, 'r') as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            parts = line.strip().lstrip('>').split('\t')
            if len(parts) < 2:
                continue

            consensus = parts[0]
            motif_id = parts[1]
            pval = '-'
            pct_targets = '-'
            pct_bg = '-'

            stats_blob = ' '.join(parts[2:])
            p_match = re.search(r'P:([0-9eE\-.]+)', stats_blob)
            t_match = re.search(r'T:[0-9.]+\(([0-9.]+%?)\)', stats_blob)
            b_match = re.search(r'B:[0-9.]+\(([0-9.]+%?)\)', stats_blob)
            if p_match:
                pval = p_match.group(1)
            if t_match:
                pct_targets = t_match.group(1)
            if b_match:
                pct_bg = b_match.group(1)

            motif_num = motif_id.split('-')[0]
            svg_candidate = homer_results_dir / f"motif{motif_num}.logo.svg"
            svg_file = str(svg_candidate) if svg_candidate.exists() else None

            motifs.append({
                'rank': len(motifs) + 1,
                'motif': motif_id,
                'consensus': consensus,
                'pval': pval,
                'pct_targets': pct_targets,
                'pct_bg': pct_bg,
                'svg_file': svg_file,
                'motif_len': len(consensus),
            })

            if len(motifs) >= top_n:
                break

    if not motifs:
        print(f"Skipping de novo PDF (no parsed motifs): {denovo_file}")
        return False

    with tempfile.TemporaryDirectory() as tmp_dir:
        for idx, m in enumerate(motifs, start=1):
            out_png = Path(tmp_dir) / f"denovo_{idx}.png"
            if _svg_to_png(m.get('svg_file'), out_png):
                m['logo_png'] = str(out_png)
            else:
                m['logo_png'] = None

        _render_logo_table_pdf(
            motifs,
            f"HOMER De Novo Motifs: {homer_dir.name} (Top {len(motifs)})",
            output_pdf,
        )

    print(f"Created motif table PDF: {output_pdf}")
    return True


def create_condition_motif_pdfs(motif_dir, top_n=5):
    """Generate known and de novo motif PDF tables per condition directory."""
    motif_root = Path(motif_dir)
    created = []

    consensus_dir = motif_root / 'consensus_peaks'
    merged_dir = motif_root / 'merged_peaks' / 'merged_peaks_motifs'

    condition_dirs = []
    if consensus_dir.exists():
        condition_dirs.extend(sorted([d for d in consensus_dir.iterdir() if d.is_dir() and d.name.endswith('_motifs')]))
    if merged_dir.exists():
        condition_dirs.append(merged_dir)

    for cond_dir in condition_dirs:
        condition = cond_dir.name.replace('_motifs', '')
        known_pdf = motif_root / f"{condition}_Known_Motifs_Table.pdf"
        denovo_pdf = motif_root / f"{condition}_DeNovo_Motifs_Table.pdf"
        if create_known_motif_pdf(cond_dir, known_pdf, top_n=top_n):
            created.append(known_pdf)
        if create_denovo_motif_pdf(cond_dir, denovo_pdf, top_n=top_n):
            created.append(denovo_pdf)

    print(f"Created {len(created)} per-condition motif PDFs")
    return created


def main():
    if len(sys.argv) < 2:
        print("Usage: create_motif_comparison_tables.py <motif_directory> [top_n]")
        print("\nExample:")
        print("  create_motif_comparison_tables.py results_dual_norm/03_peak_calling/09_homer_motifs")
        sys.exit(1)
    
    motif_dir = sys.argv[1]
    top_n = _safe_int(sys.argv[2], default=5) if len(sys.argv) > 2 else 5
    
    if not os.path.exists(motif_dir):
        print(f"Error: Directory not found: {motif_dir}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(motif_dir)
    
    print("=" * 80)
    print("HOMER MOTIF COMPARISON TABLE GENERATOR")
    print("=" * 80)
    print()
    
    # Create known motif table
    print("Processing known motifs...")
    known_table = create_known_motif_table(
        motif_dir,
        output_dir / 'Known_Motifs_Comparison_Table.tsv'
    )
    
    print()
    
    # Create de novo motif table
    print("Processing de novo motifs...")
    denovo_table = create_denovo_motif_table(
        motif_dir,
        output_dir / 'DeNovo_Motifs_Comparison_Table.tsv'
    )

    print()
    print("Creating per-condition motif PDFs...")
    create_condition_motif_pdfs(motif_dir, top_n=top_n)
    
    print()
    print("=" * 80)
    print("DONE!")
    print("=" * 80)


if __name__ == '__main__':
    main()
