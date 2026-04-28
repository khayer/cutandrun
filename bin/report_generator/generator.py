#!/usr/bin/env python3
import argparse
import json
import os
import glob
import pathlib
import pandas as pd
from jinja2 import Environment, FileSystemLoader


def find_first(glob_patterns, base):
    for pat in glob_patterns:
        matches = sorted(glob.glob(os.path.join(base, pat)))
        if matches:
            return matches
    return []


def preview_table(tsv_path, n=10):
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        return df.head(n).to_html(classes='table table-sm table-striped', index=False, na_rep='NA', escape=False)
    except Exception as e:
        return f"<pre>Error loading {os.path.basename(tsv_path)}: {e}</pre>"


def build_context(staging_dir):
    ctx = {}
    # Files of interest
    ctx['promoter_gain'] = find_first(["*.promoter_gain.tsv", "*/.promoter_gain.tsv", "**/*.promoter_gain.tsv"], staging_dir)
    ctx['promoter_loss'] = find_first(["*.promoter_loss.tsv", "**/*.promoter_loss.tsv"], staging_dir)
    ctx['deseq2'] = find_first(["*.deseq2_results.tsv", "**/*.deseq2_results.tsv"], staging_dir)
    ctx['volcano'] = find_first(["*.volcano.png", "**/*.volcano.png"], staging_dir)
    ctx['ma'] = find_first(["*.ma_plot.png", "**/*.ma_plot.png"], staging_dir)
    ctx['top_images'] = find_first(["*top*_volcano_raw_pvalue.png", "**/*top*_volcano_raw_pvalue.png", "pygt*.png", "**/pygt*.png"], staging_dir)
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

    # Convert absolute paths to relative paths from output location (we'll assume output in staging parent)
    ctx['staging_dir'] = staging_dir
    return ctx


def render_template(staging_dir, out_html):
    here = os.path.dirname(__file__)
    templates_dir = os.path.join(here, 'templates')
    env = Environment(loader=FileSystemLoader(templates_dir))
    env.filters['basename'] = os.path.basename
    tpl = env.get_template('CUT_and_RUN_report.html.j2')
    ctx = build_context(staging_dir)
    rendered = tpl.render(ctx=ctx)
    with open(out_html, 'w') as fh:
        fh.write(rendered)
    print(f"Wrote report to {out_html}")


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
