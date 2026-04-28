# CUT_and_RUN Report Generator

A Jinja2-based HTML report generator for the CUT_and_RUN pipeline that collects analysis artifacts (plots, TSVs, tables) and produces a comprehensive, publication-ready HTML report.

## Quick Start

### Local Testing

To test the report generator locally with sample data:

```bash
./run_local_report.sh
```

This will generate `tests/report_smoke/CUT_and_RUN_report.html` from the smoke test staging directory.

### Pipeline Integration

The report is generated automatically when `--with-report` is passed to the Nextflow pipeline:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --genome hg38 \
  --outdir results \
  --with-report
```

The final report will be saved to `results/report/CUT_and_RUN_report.html`.

## Files

- **generator.py** — Main report generation script. Usage: `python generator.py --staging-dir <dir> --out <output.html>`
- **templates/CUT_and_RUN_report.html.j2** — Jinja2 HTML template with sections for QC, promoter GC analysis, and annotations
- **sections.yml** — YAML mapping of section names to glob patterns for finding input artifacts
- **utils.py** — Utility functions for base64 encoding images and other helpers
- **requirements.txt** — Python dependencies (jinja2, pandas, pyyaml)
- **run_local_report.sh** — Convenience script for local testing

## Artifact Discovery

The report generator automatically scans a staging directory for pipeline outputs using glob patterns defined in `sections.yml`. Supported artifacts include:

- **Promoter GC outputs**: `*.promoter_gain.tsv`, `*.promoter_loss.tsv`, `*.deseq2_results.tsv`
- **Plots**: `*.volcano.png`, `*.ma_plot.png`, `*volcano_raw_pvalue.png`, `pygt*.png`
- **Annotations**: `*_chipseeker_annotation.tsv`, `homerResults.html`
- **QC metrics**: `*_mqc.tsv`, `frag_len_*.png`

## Customization

### Adding a New Section

1. Add a new entry to `sections.yml` with glob patterns for your artifacts
2. Extend `templates/CUT_and_RUN_report.html.j2` with a new `<div>` section using the data from the context (e.g., `ctx.your_artifact`)
3. Update `generator.py` `build_context()` function to locate and process your artifacts

### Styling

The report uses Bootstrap 5.3 CSS and DataTables for interactive tables. Modify `templates/CUT_and_RUN_report.html.j2` to adjust the layout, color scheme, or JavaScript behavior.

## Environment

The report generator runs inside a Python 3.10 container in the pipeline. Locally, it requires:

- Python 3.7+
- jinja2 ≥ 3.1.2
- pandas ≥ 1.5.0
- pyyaml ≥ 6.0

Install dependencies locally with:

```bash
pip install -r requirements.txt
```

## Troubleshooting

**No artifacts found in staging directory:**
- Ensure the staging directory path exists and contains the expected files
- Check that upstream pipeline modules are publishing their outputs with `publishDir` directives
- Run `ls -la <staging_dir>` to verify files are present

**Template rendering errors:**
- Verify all keys used in the template (e.g., `ctx.promoter_gain`) are populated by `build_context()`
- Check that Jinja2 filter usage is correct (e.g., `|basename`, `|replace`)

**Large HTML file size:**
- The report embeds full TSVs as HTML tables. For very large tables, consider publishing them separately and linking to them instead

## Future Enhancements

- Embed Plotly interactive plots for client-side zoom/hover
- Add genome browser integration (IGV snapshots)
- Support for workflow metadata (versions, command-line invocation)
- Lazy-loading of large tables
- PDF export capability
