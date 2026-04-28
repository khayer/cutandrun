#!/usr/bin/env python3
import argparse
import json
import os
import glob
import pathlib
import shutil
import pandas as pd
from jinja2 import Environment, FileSystemLoader


def find_first(glob_patterns, base):
    for pat in glob_patterns:
        # Use recursive glob so patterns with ** match nested publishDir output
        matches = sorted(glob.glob(os.path.join(base, pat), recursive=True))
        if matches:
            return matches
    return []


def preview_table(tsv_path, n=10):
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        return df.head(n).to_html(classes='table table-sm table-striped', index=False, na_rep='NA', escape=False)
    except Exception as e:
        return f"<pre>Error loading {os.path.basename(tsv_path)}: {e}</pre>"


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


def build_context(staging_dir):
    ctx = {}
    # Files of interest
    ctx['promoter_gain'] = find_first(["*.promoter_gain.tsv", "*/.promoter_gain.tsv", "**/*.promoter_gain.tsv"], staging_dir)
    ctx['promoter_loss'] = find_first(["*.promoter_loss.tsv", "**/*.promoter_loss.tsv"], staging_dir)
    ctx['deseq2'] = find_first(["*.deseq2_results.tsv", "**/*.deseq2_results.tsv"], staging_dir)
    ctx['volcano'] = find_first(["*.volcano.png", "**/*.volcano.png"], staging_dir)
    ctx['ma'] = find_first(["*.ma_plot.png", "**/*.ma_plot.png"], staging_dir)
    ctx['top_images'] = find_first(["*top*_volcano_raw_pvalue.png", "**/*top*_volcano_raw_pvalue.png", "pygt*.png", "**/pygt*.png"], staging_dir)
    # Find pygenometracks directories grouped by condition
    pygt_base = sorted(glob.glob(os.path.join(staging_dir, '**', 'pygenometracks_top10'), recursive=True))
    ctx['pygenometracks_dirs'] = {}
    if pygt_base:
        for pygt_dir in pygt_base:
            for condition_dir in glob.glob(os.path.join(pygt_dir, '*')):
                if os.path.isdir(condition_dir):
                    condition_name = os.path.basename(condition_dir)
                    images = sorted(glob.glob(os.path.join(condition_dir, '*.png')))
                    if images:
                        ctx['pygenometracks_dirs'][condition_name] = images
    ctx['chipseeker'] = find_first(["*_chipseeker_annotation.tsv", "**/*_chipseeker_annotation.tsv"], staging_dir)
    ctx['homer_html'] = find_first(["**/homerResults.html", "homerResults.html"], staging_dir)
    ctx['frip'] = find_first(["*_mqc.tsv", "**/*_mqc.tsv"], staging_dir)

    # Preview HTMLs
    ctx['previews'] = {}
    if ctx['promoter_gain']:
        ctx['previews']['promoter_gain'] = preview_table(ctx['promoter_gain'][0])
    if ctx['promoter_loss']:
        ctx['previews']['promoter_loss'] = preview_table(ctx['promoter_loss'][0])
    if ctx['deseq2']:
        ctx['previews']['deseq2'] = preview_table(ctx['deseq2'][0])
    if ctx['chipseeker']:
        ctx['previews']['chipseeker'] = preview_table(ctx['chipseeker'][0])
    
    # FRiP scores table
    if ctx['frip']:
        ctx['previews']['frip'] = parse_frip_scores(ctx['frip'])

    # Convert absolute paths to relative paths from output location (we'll assume output in staging parent)
    ctx['staging_dir'] = staging_dir
    return ctx


def render_template(staging_dir, out_html):
    here = os.path.dirname(__file__)
    templates_dir = os.path.join(here, 'templates')
    env = Environment(loader=FileSystemLoader(templates_dir))
    env.filters['basename'] = os.path.basename
    env.filters['dirname'] = os.path.dirname
    tpl = env.get_template('CUT_and_RUN_report.html.j2')
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
        dest_name = os.path.basename(src_path)
        dest_path = os.path.join(artifacts_dir, dest_name)
        if not os.path.exists(dest_path):
            shutil.copy2(src_path, dest_path)
        return f"artifacts/{dest_name}"
    
    # Convert artifact lists to links within report folder
    ctx['promoter_gain_rel'] = [copy_artifact(p) for p in (ctx.get('promoter_gain') or [])]
    ctx['promoter_loss_rel'] = [copy_artifact(p) for p in (ctx.get('promoter_loss') or [])]
    ctx['deseq2_rel'] = [copy_artifact(p) for p in (ctx.get('deseq2') or [])]
    ctx['volcano_rel'] = [copy_artifact(p) for p in (ctx.get('volcano') or [])]
    ctx['ma_rel'] = [copy_artifact(p) for p in (ctx.get('ma') or [])]
    ctx['chipseeker_rel'] = [copy_artifact(p) for p in (ctx.get('chipseeker') or [])]
    ctx['homer_html_rel'] = [copy_artifact(p) for p in (ctx.get('homer_html') or [])]
    
    # Copy pygenometracks images organized by condition
    ctx['pygenometracks_dirs_rel'] = {}
    for condition, img_list in ctx.get('pygenometracks_dirs', {}).items():
        ctx['pygenometracks_dirs_rel'][condition] = [copy_artifact(img) for img in img_list]
    
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
    args = parse_args()
    staging = os.path.abspath(args.staging_dir)
    out = os.path.abspath(args.out)
    if not os.path.isdir(staging):
        raise SystemExit(f"Staging dir not found: {staging}")
    render_template(staging, out)
