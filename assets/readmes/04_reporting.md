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
Enabled with `--run_pygenometracks_top10 true`. For each non-control group, the top `--pygt_top_n` consensus peaks (ranked by consensus support score) are rendered as high-resolution track plots using [pyGenomeTracks](https://pygenometracks.readthedocs.io/). Each panel shows:
- All non-control replicates for the focal group
- One IgG/control track (grey)
- One representative replicate from each other group (for comparison)

Output per group (`<GROUP>/`):
| File | Description |
|------|-------------|
| `top_regions.tsv` | Selected top regions with nearest-gene annotation and plot windows |
| `top_peaks.bed` | BED file of selected top peaks |
| `tracks.ini` | pyGenomeTracks configuration used for rendering |
| `<GROUP>_top##_*.png` / `.pdf` | Rendered locus snapshots |

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
