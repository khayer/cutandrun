# 04_reporting — Quality Control Reports and Visualisations

This directory contains all quality-control visualisations and summary reports produced by the pipeline.

## Subdirectories

| Path | Contents |
|------|----------|
| `multiqc/` | **MultiQC interactive HTML report** — aggregates all QC metrics in one page |
| `deeptools_qc/` | deepTools fragment-level QC: PCA, fingerprint, and correlation matrices |
| `deeptools_heatmaps/gene/` | deepTools heatmap of fragment coverage over gene bodies |
| `deeptools_heatmaps/gene_all/` | Gene-body heatmap including all samples |
| `deeptools_heatmaps/peaks/` | deepTools heatmap of fragment coverage centred on called peaks |
| `igv/` | IGV session XML file linking normalised bigWig tracks and peaks |
| `igv_downsampled/` | IGV session using downsampled (equal-depth) bigWig tracks |
| `pygenometracks_top10/` | pyGenomeTracks locus snapshots for the top consensus peaks per group (when `--run_pygenometracks_top10 true`) |
| `preseq/` | Preseq library complexity curves |
| `consensus_upset_plots/` | Upset plots showing peak-set overlap across groups |
| `peaksignalprofiler/` | PeakSignalProfiler outputs (when `--run_peak_signal_profiler true`) |

## How These Files Are Derived

### MultiQC (`multiqc/`)
[MultiQC](https://multiqc.info/) aggregates log files and metrics from every tool run in the pipeline into a single interactive HTML report (`multiqc_report.html`). It includes:
- FastQC sequence quality and adapter content (pre- and post-trim)
- Bowtie2 alignment rates (target and spike-in)
- Picard duplication metrics
- Preseq library complexity
- deepTools PCA, fingerprint, and correlation
- Peak counts and FRiP scores
- Fragment length distribution histogram

**Start here when evaluating a new run.**

### deepTools QC (`deeptools_qc/`)
[deepTools](https://deeptools.readthedocs.io/) bins the genome into windows (default: 500 bp; configurable with `--dt_qc_bam_binsize`) and counts fragments per window across all samples. Three analyses are performed:

- **PCA** — Principal Component Analysis on the fragment-count matrix. Samples from the same experimental group should cluster together. Outliers or unexpected groupings indicate batch effects or sample swap.
- **Fingerprint** — Plots the cumulative fraction of reads in increasingly enriched bins. A diagonal line (uniform coverage) indicates no enrichment; a strongly curved line indicates a well-enriched sample.
- **Correlation** — Pairwise Spearman or Pearson correlation of fragment counts between all samples. Biological replicates should show high correlation (≥ 0.9 for histone marks).

### deepTools Heatmaps (`deeptools_heatmaps/`)
Coverage heatmaps are generated using `deepTools computeMatrix` and `plotHeatmap`:
- **Gene heatmaps** — Fragment coverage from `--dt_heatmap_gene_beforelen` kb before TSS to `--dt_heatmap_gene_afterlen` kb after TES, scaled to gene body length.
- **Peak heatmaps** — Fragment coverage centred on the called peaks, flanked by `--dt_heatmap_peak_beforelen` / `--dt_heatmap_peak_afterlen` bp.

### IGV Session (`igv/`, `igv_downsampled/`)
An [IGV](https://igv.org/) XML session file is generated that pre-loads:
- Normalised bigWig tracks for all samples
- Per-sample peak BED files
- The target genome FASTA and GTF annotation

To use: open IGV, go to **File → Open Session** and select `igv_session.xml`.

### pyGenomeTracks Top Peaks (`pygenometracks_top10/`)
Enabled with `--run_pygenometracks_top10 true`. For each non-control experimental group, the pipeline selects the top N consensus peaks (ranked by consensus support score — the number of replicates in which the peak was called) and renders them as high-resolution genomic track images using [pyGenomeTracks](https://pygenometracks.readthedocs.io/).

**How top peaks are selected:**
1. Consensus peaks for the group are loaded from the `05_consensus_peaks/` BED files.
2. Peaks are ranked by their score column (consensus support from the BED file, column 10 if present, otherwise column 5).
3. The top `--pygt_top_n` peaks (default: 10) are selected.
4. For each peak, the midpoint is identified and a window of ± `--pygt_peak_flank` bp (default: 1,000 bp) is defined as the plot region.
5. The nearest gene to each peak midpoint is identified from the GTF and included in the plot title.

**Track composition per output image:**
- All non-control replicates for the focal group (one bigWig track each)
- One representative IgG/control track (rendered in grey)
- One representative replicate from each other non-control group (for cross-group comparison)
- A gene annotation track from the GTF
- A consensus peaks BED track

**Output per group** (one subdirectory `<GROUP>/` per non-control group):

| File | Description |
|------|-------------|
| `top_regions.tsv` | Tab-delimited table of selected regions: rank, coordinates, plot window, nearest gene name, and distance to gene |
| `top_peaks.bed` | BED file of the selected top peaks |
| `tracks.ini` | pyGenomeTracks configuration file used for rendering |
| `<GROUP>_top##_<chrom>_<start>_<end>.png` / `.pdf` | Rendered locus snapshot for each selected peak |

**Relevant parameters:**
- `--run_pygenometracks_top10` (default: false) — enable top-peak snapshots
- `--pygt_top_n` (default: 10) — number of peaks to render per group
- `--pygt_peak_flank` (default: 1000) — flanking region in bp around each peak midpoint
- `--pygt_output_format` (default: `png`) — output format: `png`, `pdf`, or `both`

### Library Complexity (`preseq/`)
[Preseq](http://smithlabresearch.org/software/preseq/) estimates how many additional unique reads would be gained with more sequencing. The output curves are included in MultiQC. A rapidly saturating curve means the library is near saturation and further sequencing will yield diminishing returns.

### Upset Plots (`consensus_upset_plots/`)
[Upset plots](https://upset.app/) visualise intersections between consensus peak sets from different experimental groups. They are more informative than Venn diagrams when ≥ 3 groups are compared.

### PeakSignalProfiler (`peaksignalprofiler/`)
Enabled with `--run_peak_signal_profiler true`. PeakSignalProfiler (PSP) calculates enrichment metrics for each sample relative to the consensus peaks and produces summary plots. See the PSP documentation for detailed interpretation.

## Recommended Interpretation Checklist

- [ ] Open `multiqc/multiqc_report.html` and check all summary statistics
- [ ] Verify alignment rates are ≥ 80% for target genome
- [ ] Check spike-in alignment levels are consistent across samples
- [ ] Confirm duplication rates are similar between biological replicates
- [ ] Inspect PCA — replicates should cluster, groups should separate
- [ ] Check FRiP scores ≥ 0.3
- [ ] Open `igv/igv_session.xml` in IGV and visually inspect peaks at known loci
- [ ] Review heatmaps for expected enrichment patterns
- [ ] Inspect `pygenometracks_top10/<GROUP>/` snapshots to visually confirm top peaks look biologically credible at individual loci (requires `--run_pygenometracks_top10 true`)
