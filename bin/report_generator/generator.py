#!/usr/bin/env python3
import argparse
import csv
import json
import os
import glob
import fnmatch
import pathlib
import shutil
import pandas as pd
from jinja2 import Environment, FileSystemLoader


def find_first(glob_patterns, base):
    # Directories to skip during recursive search
    skip_dirs = {'work', '.git', 'logs', 'tmp', '__pycache__', '.nextflow'}
    
    for pat in glob_patterns:
        if '**' not in pat:
            # Non-recursive pattern - use glob as-is
            matches = sorted(glob.glob(os.path.join(base, pat), recursive=False))
            if matches:
                return matches
        else:
            # Recursive pattern - use walk to have better control
            matches = []
            for root, dirs, files in os.walk(base):
                # Filter out skip_dirs in-place to prevent traversal
                dirs[:] = [d for d in dirs if d not in skip_dirs]
                
                # Convert pattern to filename pattern for matching
                pat_filename = pat.replace('**/', '').lstrip('/')
                for file in files:
                    if fnmatch.fnmatch(file, pat_filename):
                        matches.append(os.path.join(root, file))
            
            if matches:
                return sorted(matches)
    
    return []


def preview_table(path, n=10):
    try:
        # Auto-detect separator for CSV vs TSV
        if path.endswith('.csv'):
            df = pd.read_csv(path, dtype=str)
        else:
            df = pd.read_csv(path, sep="\t", dtype=str)
        # Show only first 10 columns if there are many
        if len(df.columns) > 10:
            df = df.iloc[:, :10]
        return df.head(n).to_html(classes='table table-sm table-striped', index=False, na_rep='NA', escape=False)
    except Exception as e:
        return f"<pre>Error loading {os.path.basename(path)}: {e}</pre>"


def parse_frip_scores(frip_files):
    """Parse FRiP TSV files and extract scores into a table."""
    if not frip_files:
        return None
    data = []
    for fpath in frip_files:
        try:
            sample_name = os.path.basename(fpath).replace('_mqc.tsv', '')
            with open(fpath, 'r') as f:
                for line in f:
                    if line.startswith('Peak FRiP Score'):
                        parts = line.strip().split()
                        score = parts[-1]
                        data.append({'Sample': sample_name, 'FRiP Score': score})
                        break
        except Exception as e:
            pass
    if data:
        df = pd.DataFrame(data)
        return df.to_html(classes='table table-sm table-striped', index=False, na_rep='NA', escape=False)
    return None


def parse_gc_test(gc_test_path):
    """Parse gc_test.tsv file and return as HTML table."""
    try:
        df = pd.read_csv(gc_test_path, sep="\t")
        # Format numeric columns for display
        if 'p_value' in df.columns:
            df['p_value'] = df['p_value'].apply(lambda x: f"{x:.2e}" if pd.notna(x) else x)
        if 'statistic' in df.columns:
            df['statistic'] = df['statistic'].apply(lambda x: f"{x:.2f}" if pd.notna(x) else x)
        for col in ['median_gc_loss', 'median_gc_gain', 'median_difference']:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else x)
        return df.to_html(classes='table table-sm table-striped', index=False, na_rep='NA', escape=False)
    except Exception as e:
        return f"<pre>Error loading gc_test.tsv: {e}</pre>"


def build_context(staging_dir):
    ctx = {}
    def is_report_artifact(path):
        return f"{os.path.sep}report{os.path.sep}artifacts{os.path.sep}" in path.replace('/', os.path.sep)

    # Files of interest
    ctx['promoter_gain'] = find_first(["*.promoter_gain.tsv", "*/.promoter_gain.tsv", "**/*.promoter_gain.tsv"], staging_dir)
    ctx['promoter_loss'] = find_first(["*.promoter_loss.tsv", "**/*.promoter_loss.tsv"], staging_dir)
    ctx['deseq2'] = find_first(["*.deseq2_results.tsv", "**/*.deseq2_results.tsv"], staging_dir)
    
    # DESeq2 peak analysis results (CSV format)
    deseq2_results_files = sorted(glob.glob(os.path.join(staging_dir, '**', '*_results.csv'), recursive=True))
    deseq2_results_files = [f for f in deseq2_results_files if not is_report_artifact(f)]
    ctx['deseq2_results'] = deseq2_results_files[0] if deseq2_results_files else None  # Take first if multiple
    
    deseq2_sig_files = sorted(glob.glob(os.path.join(staging_dir, '**', '*_significant.csv'), recursive=True))
    deseq2_sig_files = [f for f in deseq2_sig_files if not is_report_artifact(f)]
    ctx['deseq2_significant'] = deseq2_sig_files[0] if deseq2_sig_files else None  # Take first if multiple
    
    # DESeq2 plots - PNG preferred, fallback to PDF
    deseq2_volcano_files = sorted(glob.glob(os.path.join(staging_dir, '**', '*_volcano_plot.png'), recursive=True))
    deseq2_volcano_files = [f for f in deseq2_volcano_files if not is_report_artifact(f)]
    if not deseq2_volcano_files:
        deseq2_volcano_files = sorted(glob.glob(os.path.join(staging_dir, '**', '*_volcano_plot.pdf'), recursive=True))
        deseq2_volcano_files = [f for f in deseq2_volcano_files if not is_report_artifact(f)]
    ctx['deseq2_volcano'] = deseq2_volcano_files[0] if deseq2_volcano_files else None
    ctx['deseq2_volcano_is_png'] = ctx['deseq2_volcano'].endswith('.png') if ctx['deseq2_volcano'] else False
    
    deseq2_ma_files = sorted(glob.glob(os.path.join(staging_dir, '**', '*_ma_plot.png'), recursive=True))
    deseq2_ma_files = [f for f in deseq2_ma_files if not is_report_artifact(f)]
    if not deseq2_ma_files:
        deseq2_ma_files = sorted(glob.glob(os.path.join(staging_dir, '**', '*_ma_plot.pdf'), recursive=True))
        deseq2_ma_files = [f for f in deseq2_ma_files if not is_report_artifact(f)]
    ctx['deseq2_ma'] = deseq2_ma_files[0] if deseq2_ma_files else None
    ctx['deseq2_ma_is_png'] = ctx['deseq2_ma'].endswith('.png') if ctx['deseq2_ma'] else False
    
    ctx['volcano'] = find_first(["*.volcano.png", "**/*.volcano.png"], staging_dir)
    ctx['ma'] = find_first(["*.ma_plot.png", "**/*.ma_plot.png"], staging_dir)
    ctx['top_images'] = find_first(["*top*_volcano_raw_pvalue.png", "**/*top*_volcano_raw_pvalue.png", "pygt*.png", "**/pygt*.png"], staging_dir)
    # Find pygenometracks directories grouped by condition and feature
    pygt_base = [path for path in sorted(glob.glob(os.path.join(staging_dir, '**', 'pygenometracks_top10'), recursive=True)) if not is_report_artifact(path)]
    ctx['pygenometracks'] = {}
    if pygt_base:
        for pygt_dir in pygt_base:
            for condition_dir in sorted(glob.glob(os.path.join(pygt_dir, '*'))):
                if os.path.isdir(condition_dir):
                    condition_name = os.path.basename(condition_dir)
                    feature_groups = {}
                    for image_path in sorted(glob.glob(os.path.join(condition_dir, '*.png'))):
                        basename = os.path.basename(image_path)
                        prefix = f"{condition_name}_"
                        remainder = basename[len(prefix):] if basename.startswith(prefix) else basename
                        parts = remainder.split('_')
                        feature = 'Other'
                        if len(parts) >= 1:
                            feature_token = parts[0].lower()
                            if feature_token == 'tes':
                                feature = 'TES'
                            else:
                                feature = feature_token.capitalize()
                        feature_groups.setdefault(feature, []).append(image_path)
                    if feature_groups:
                        ctx['pygenometracks'][condition_name] = feature_groups
    ctx['chipseeker'] = find_first(["*_chipseeker_annotation.tsv", "**/*_chipseeker_annotation.tsv"], staging_dir)
    ctx['homer_html'] = find_first(["**/homerResults.html", "homerResults.html"], staging_dir)
    # PNG-only mode: skip Homer table PDFs
    ctx['homer_tables'] = []
    
    # Peak Feature Annotation GC Content Plots
    gc_plot_patterns = [
        "homer_peak_annotation.gc_by_sample.png",
        "homer_peak_annotation.gc_by_sample_violin.png",
        "homer_peak_annotation.gc_by_condition.png",
        "homer_peak_annotation.gc_by_condition_violin.png",
    ]
    ctx['gc_plots'] = {}
    for pattern in gc_plot_patterns:
        files = find_first([f"*/{pattern}", f"**/{pattern}"], staging_dir)
        if files:
            ctx['gc_plots'][pattern.replace('homer_peak_annotation.', '').replace('.png', '')] = files[0]
    
    # ChIPseeker Coverage Comparison Plots
    coverage_patterns = [
        "chipseeker_comparison_coverage.png",
        "chipseeker_comparison_coverage_condition_average.png",
    ]
    ctx['coverage_plots'] = {}
    for pattern in coverage_patterns:
        files = find_first([f"*/{pattern}", f"**/{pattern}"], staging_dir)
        if files:
            label = pattern.replace('chipseeker_comparison_coverage', '').replace('_condition_average', ' (by condition)').replace('.png', '').lstrip('_') or 'all samples'
            ctx['coverage_plots'][label] = files[0]
    
    # Per-Replicate ChIPseeker Plots - Group by condition and type
    ctx['per_replicate_chipseeker'] = {}
    all_chipseeker_pngs = []
    for pattern in ["**/*_R1_chipseeker*.png", "**/*_R2_chipseeker*.png", "**/*_R3_chipseeker*.png"]:
        all_chipseeker_pngs.extend(glob.glob(os.path.join(staging_dir, pattern), recursive=True))
    all_chipseeker_pngs = [p for p in all_chipseeker_pngs if not is_report_artifact(p)]
    
    for png_path in sorted(all_chipseeker_pngs):
        basename = os.path.basename(png_path)
        # Extract condition and replicate: e.g., "33K_Ad5_R1_chipseeker_annotation_pie.png"
        # becomes condition="33K_Ad5", replicate="R1", plot_type="annotation_pie"
        parts = basename.replace('.png', '').split('_')
        if '_R' in basename:
            r_idx = -1
            for i, p in enumerate(parts):
                if p.startswith('R') and p[1:].isdigit():
                    r_idx = i
                    break
            if r_idx >= 0:
                condition = '_'.join(parts[:r_idx])
                replicate = parts[r_idx]
                plot_type = '_'.join(parts[r_idx+2:])  # Skip 'chipseeker' part
                
                if condition not in ctx['per_replicate_chipseeker']:
                    ctx['per_replicate_chipseeker'][condition] = {}
                if replicate not in ctx['per_replicate_chipseeker'][condition]:
                    ctx['per_replicate_chipseeker'][condition][replicate] = {}
                ctx['per_replicate_chipseeker'][condition][replicate][plot_type] = png_path
    
    # Alternative ChIPseeker Comparison Plots (non-standard variants)
    ctx['alternative_chipseeker_plots'] = {}
    alternative_patterns = [
        'chipseeker_comparison_annotation_bar_slim_version',
        'chipseeker_comparison_annotation_bar_condition',
        'chipseeker_comparison_annotation_bar_condition_slim_version',
        'chipseeker_condition_.*_pie',
    ]
    for pattern in alternative_patterns:
        files = sorted(glob.glob(os.path.join(staging_dir, '**', f'{pattern}.png'), recursive=True))
        files = [f for f in files if not is_report_artifact(f)]
        if files:
            label = pattern.replace('chipseeker_', '').replace('_', ' ').title()
            for f in files:
                ctx['alternative_chipseeker_plots'][f"{label} ({os.path.basename(f)})"] = f
    ctx['diffbind_plots'] = {}
    ctx['diffbind_plots_pval'] = {}  # Raw p-value version plots
    ctx['diffbind_data'] = {}  # Store gain/loss/gc_test files per contrast
    ctx['diffbind_data_pval'] = {}  # Raw p-value version data
    diffbind_roots = [path for path in sorted(glob.glob(os.path.join(staging_dir, '**', '11_differential_binding', '*'))) if not is_report_artifact(path)]
    plot_priority = [
        'top1000_promoter_gc_boxplot.png',
        'top1000_promoter_gc_violin.png',
        'top1000_volcano_raw_pvalue.png',
        'promoter_gc_boxplot.png',
        'promoter_gc_violin.png',
        'volcano.png',
        'ma_plot.png',
    ]
    plot_priority_pval = [
        'promoter_gc_boxplot_raw_pval.png',
        'promoter_gc_violin_raw_pval.png',
    ]
    for contrast_dir in diffbind_roots:
        if not os.path.isdir(contrast_dir):
            continue
        if os.path.join('report', 'artifacts') in contrast_dir.replace(os.path.sep, '/'):
            continue
        contrast_name = os.path.basename(contrast_dir)
        plot_map = {}
        for image_path in glob.glob(os.path.join(contrast_dir, '*.png')):
            plot_map[os.path.basename(image_path)] = image_path
        ordered_plots = []
        used = set()
        for priority_suffix in plot_priority:
            for filename, image_path in plot_map.items():
                if filename.endswith(priority_suffix) and filename not in used:
                    ordered_plots.append(image_path)
                    used.add(filename)
        for filename, image_path in sorted(plot_map.items()):
            if filename not in used:
                ordered_plots.append(image_path)
        if ordered_plots:
            ctx['diffbind_plots'][contrast_name] = ordered_plots
        
        # Find top 1000 gain/loss/gc_test files for this contrast
        gain_file = None
        loss_file = None
        gc_test_file = None
        for file_path in glob.glob(os.path.join(contrast_dir, '*.top1000_gain.tsv')):
            gain_file = file_path
            break
        for file_path in glob.glob(os.path.join(contrast_dir, '*.top1000_loss.tsv')):
            loss_file = file_path
            break
        for file_path in glob.glob(os.path.join(contrast_dir, '*.top1000_promoter_gc_test.tsv')):
            gc_test_file = file_path
            break
        if gain_file or loss_file or gc_test_file:
            ctx['diffbind_data'][contrast_name] = {
                'gain': gain_file,
                'loss': loss_file,
                'gc_test': gc_test_file,
            }
        
        # Find raw p-value gc_test files for this contrast
        gc_test_pval_file = None
        for file_path in glob.glob(os.path.join(contrast_dir, '*.promoter_gc_test_raw_pval.tsv')):
            gc_test_pval_file = file_path
            break
        if gc_test_pval_file:
            ctx['diffbind_data_pval'][contrast_name] = {
                'gc_test': gc_test_pval_file,
            }
        
        # Collect raw p-value plots
        pval_plots = []
        for priority_suffix in plot_priority_pval:
            for filename, image_path in plot_map.items():
                if filename.endswith(priority_suffix):
                    pval_plots.append(image_path)
        if pval_plots:
            ctx['diffbind_plots_pval'][contrast_name] = pval_plots
    # Prefer PNG for faster loading, fall back to PDF/JPG
    chp_plots = []
    for pattern in ["**/*chipseeker*png", "**/*chipseeker*jpg", "**/*chipseeker*pdf"]:
        chp_plots.extend(path for path in glob.glob(os.path.join(staging_dir, pattern), recursive=True) if not is_report_artifact(path))
    ctx['chipseeker_plots'] = sorted(set(chp_plots))
    ctx['frip'] = find_first(["*_mqc.tsv", "**/*_mqc.tsv"], staging_dir)

    # Preview HTMLs
    ctx['previews'] = {}
    if ctx['promoter_gain']:
        ctx['previews']['promoter_gain'] = preview_table(ctx['promoter_gain'][0])
    if ctx['promoter_loss']:
        ctx['previews']['promoter_loss'] = preview_table(ctx['promoter_loss'][0])
    if ctx['deseq2']:
        ctx['previews']['deseq2'] = preview_table(ctx['deseq2'][0])
    if ctx['deseq2_significant']:
        ctx['previews']['deseq2_significant'] = preview_table(ctx['deseq2_significant'])
    if ctx['chipseeker']:
        ctx['previews']['chipseeker'] = preview_table(ctx['chipseeker'][0])
    
    # FRiP scores table
    if ctx['frip']:
        ctx['previews']['frip'] = parse_frip_scores(ctx['frip'])

    # Convert absolute paths to relative paths from output location (we'll assume output in staging parent)
    ctx['staging_dir'] = staging_dir
    return ctx


def render_template(staging_dir, out_html):
    print(f"[DEBUG] render_template called with staging_dir={staging_dir}", flush=True)
    print(f"[DEBUG] Checking staging directory exists...", flush=True)
    if not os.path.exists(staging_dir):
        raise ValueError(f"Staging directory does not exist: {staging_dir}")
    print(f"[DEBUG] Staging directory OK", flush=True)
    
    here = os.path.dirname(__file__)
    templates_dir = os.path.join(here, 'templates')
    print(f"[DEBUG] Loading Jinja2 templates from: {templates_dir}", flush=True)
    
    env = Environment(loader=FileSystemLoader(templates_dir))
    env.filters['basename'] = os.path.basename
    env.filters['dirname'] = os.path.dirname
    tpl = env.get_template('CUT_and_RUN_report.html.j2')
    
    print(f"[DEBUG] Building context from staging_dir...", flush=True)
    ctx = build_context(staging_dir)
    
    # Setup output directory and artifacts folder
    out_dir = os.path.dirname(out_html)
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    
    artifacts_dir = os.path.join(out_dir, 'artifacts')
    if not os.path.exists(artifacts_dir):
        os.makedirs(artifacts_dir, exist_ok=True)
    
    # Copy artifacts into report folder for portability
    def copy_artifact(src_path):
        """Copy artifact to artifacts/ folder and return relative link."""
        if not src_path or not os.path.exists(src_path):
            return None
        rel_src = os.path.relpath(src_path, start=staging_dir)
        rel_src = rel_src.replace(os.path.sep, '/')
        dest_path = os.path.join(artifacts_dir, rel_src)
        dest_dir = os.path.dirname(dest_path)
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir, exist_ok=True)
        if not os.path.exists(dest_path):
            shutil.copy2(src_path, dest_path)
        return f"artifacts/{rel_src}"
    
    # Convert artifact lists to links within report folder
    ctx['promoter_gain_rel'] = [copy_artifact(p) for p in (ctx.get('promoter_gain') or [])]
    ctx['promoter_loss_rel'] = [copy_artifact(p) for p in (ctx.get('promoter_loss') or [])]
    ctx['deseq2_rel'] = [copy_artifact(p) for p in (ctx.get('deseq2') or [])]
    ctx['deseq2_results_rel'] = copy_artifact(ctx.get('deseq2_results')) if ctx.get('deseq2_results') else None
    ctx['deseq2_significant_rel'] = copy_artifact(ctx.get('deseq2_significant')) if ctx.get('deseq2_significant') else None
    
    # PNG-only mode for all plots
    ctx['deseq2_volcano_rel'] = copy_artifact(ctx.get('deseq2_volcano')) if ctx.get('deseq2_volcano') else None
    ctx['deseq2_ma_rel'] = copy_artifact(ctx.get('deseq2_ma')) if ctx.get('deseq2_ma') else None
    # Note: is_png flags already set in build_context, preserve them here
    ctx['volcano_rel'] = [copy_artifact(p) for p in (ctx.get('volcano') or [])]
    ctx['ma_rel'] = [copy_artifact(p) for p in (ctx.get('ma') or [])]
    ctx['chipseeker_rel'] = [copy_artifact(p) for p in (ctx.get('chipseeker') or [])]
    ctx['homer_reports'] = []
    for p in (ctx.get('homer_html') or []):
        label = os.path.basename(os.path.dirname(p))
        ctx['homer_reports'].append({'label': label, 'link': copy_artifact(p)})
    # PNG-only mode: skip Homer table PDF artifacts
    ctx['homer_tables_rel'] = []
    
    # GC Content Plots
    ctx['gc_plots_rel'] = {}
    for label, path in ctx.get('gc_plots', {}).items():
        ctx['gc_plots_rel'][label] = copy_artifact(path)
    
    # Coverage Comparison Plots
    ctx['coverage_plots_rel'] = {}
    for label, path in ctx.get('coverage_plots', {}).items():
        ctx['coverage_plots_rel'][label] = copy_artifact(path)
    
    # Per-Replicate ChIPseeker Plots
    ctx['per_replicate_chipseeker_rel'] = {}
    for condition, replicates in ctx.get('per_replicate_chipseeker', {}).items():
        ctx['per_replicate_chipseeker_rel'][condition] = {}
        for replicate, plots in replicates.items():
            ctx['per_replicate_chipseeker_rel'][condition][replicate] = {}
            for plot_type, path in plots.items():
                ctx['per_replicate_chipseeker_rel'][condition][replicate][plot_type] = copy_artifact(path)
    
    # Alternative ChIPseeker Comparison Plots
    ctx['alternative_chipseeker_plots_rel'] = {}
    for label, path in ctx.get('alternative_chipseeker_plots', {}).items():
        ctx['alternative_chipseeker_plots_rel'][label] = copy_artifact(path)
    ctx['chipseeker_comparison_plots_rel'] = [
        copy_artifact(p)
        for p in (ctx.get('chipseeker_plots') or [])
        if os.path.basename(p).startswith('chipseeker_comparison') and os.path.basename(p).endswith('.png')
    ]
    ctx['diffbind_plots_rel'] = {}
    for contrast, plots in ctx.get('diffbind_plots', {}).items():
        ctx['diffbind_plots_rel'][contrast] = [copy_artifact(plot) for plot in plots]
    ctx['diffbind_plots_pval_rel'] = {}
    for contrast, plots in ctx.get('diffbind_plots_pval', {}).items():
        ctx['diffbind_plots_pval_rel'][contrast] = [copy_artifact(plot) for plot in plots]
    ctx['diffbind_data_rel'] = {}
    for contrast, data in ctx.get('diffbind_data', {}).items():
        gc_test_html = None
        if data.get('gc_test'):
            gc_test_html = parse_gc_test(data['gc_test'])
        ctx['diffbind_data_rel'][contrast] = {
            'gain': copy_artifact(data['gain']) if data.get('gain') else None,
            'loss': copy_artifact(data['loss']) if data.get('loss') else None,
            'gc_test': copy_artifact(data['gc_test']) if data.get('gc_test') else None,
            'gc_test_html': gc_test_html,
        }
    ctx['diffbind_data_pval_rel'] = {}
    for contrast, data in ctx.get('diffbind_data_pval', {}).items():
        gc_test_pval_html = None
        if data.get('gc_test'):
            gc_test_pval_html = parse_gc_test(data['gc_test'])
        ctx['diffbind_data_pval_rel'][contrast] = {
            'gc_test': copy_artifact(data['gc_test']) if data.get('gc_test') else None,
            'gc_test_html': gc_test_pval_html,
        }
    
    # Copy pygenometracks images organized by condition and feature
    ctx['pygenometracks_rel'] = {}
    for condition, feats in ctx.get('pygenometracks', {}).items():
        ctx['pygenometracks_rel'][condition] = {}
        for feature, img_list in feats.items():
            items = []
            for img in img_list:
                # PNG-only mode: skip PDF alternatives
                items.append({
                    'png': copy_artifact(img),
                    'pdf': None,  # Removed PDF alternatives for PNG-only report
                })
            ctx['pygenometracks_rel'][condition][feature] = items
    
    rendered = tpl.render(ctx=ctx)
    with open(out_html, 'w') as fh:
        fh.write(rendered)
    
    # Copy static assets (CSS, JS) to output directory for portability
    static_src = os.path.join(here, 'static')
    static_dest = os.path.join(out_dir, 'static')
    if os.path.isdir(static_src):
        if os.path.exists(static_dest):
            shutil.rmtree(static_dest)
        shutil.copytree(static_src, static_dest)
    
    print(f"Wrote report to {out_html}")
    print(f"Artifacts copied to {artifacts_dir}")


def parse_args():
    p = argparse.ArgumentParser(description='Generate CUT_and_RUN_report.html from staging directory')
    p.add_argument('--staging-dir', required=True)
    p.add_argument('--out', required=True)
    return p.parse_args()


if __name__ == '__main__':
    try:
        args = parse_args()
        staging = os.path.abspath(args.staging_dir)
        out = os.path.abspath(args.out)
        if not os.path.isdir(staging):
            raise SystemExit(f"Staging dir not found: {staging}")
        print(f"[INFO] Generating report from: {staging}", flush=True)
        render_template(staging, out)
        print(f"[SUCCESS] Report written to: {out}", flush=True)
    except Exception as e:
        import traceback
        print(f"[ERROR] Report generation failed: {e}", flush=True)
        print(traceback.format_exc(), flush=True)
        raise SystemExit(f"Report generation failed: {e}")
