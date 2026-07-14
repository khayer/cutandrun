# 03_peak_calling — Coverage Tracks, Peak Calling, and Downstream Analysis

This directory contains all outputs related to generating coverage tracks, calling enrichment peaks, and performing downstream analyses such as motif discovery and peak annotation.

## Subdirectories

| Path | Contents |
|------|----------|
| `01_bam_to_bedgraph/` | Per-sample bedgraph coverage files (intermediate; used to generate bigWig and call peaks) |
| `03_bed_to_bigwig/` | **Normalised** bigWig coverage tracks — primary files for genome-browser visualisation |
| `03_bed_to_bigwig_downsampled/` | Downsampled bigWig tracks (equal read depth across samples; for visual comparison) |
| `04_called_peaks/seacr/` | Per-sample peaks called by SEACR |
| `04_called_peaks/macs2/` | Per-sample peaks called by MACS2 |
| `05_consensus_peaks/` | Per-group consensus peaks (union of replicate peaks) |
| `06_fragments_from_bams/` | Fragment BED files extracted from BAMs (used for FRiP score calculation) |
| `07_bigwig_minus_igg/` | BigWig tracks with IgG signal subtracted (for background-corrected visualisation) |
| `07_peak_qc/frip/` | Fraction of Reads in Peaks (FRiP) score tables |
| `08_merged_peaks_table/` | Merged peak table with read counts across all samples |
| `09_homer_motifs/` | De-novo and known motif discovery results (only when `--run_homer_motifs true`) |
| `10_peak_feature_annotation/` | Genomic feature annotation of peaks (HOMER summary + ChIPseeker visualisations; requires `--run_homer_peak_annotation true`) |

## How These Files Are Derived

### Coverage Tracks: bedgraph → bigWig

1. **bedgraph** (`01_bam_to_bedgraph/`) — BAM files are converted to bedgraph format using `bedtools genomecov`, computing the number of fragments (or reads) overlapping each base. Spike-in normalisation is applied at this stage: the bedgraph values are scaled by `C / (spike-in_aligned_reads)`, where `C` is a constant (`--normalisation_c`, default: 10,000), so that samples with different spike-in fractions are brought to the same scale.

2. **bigWig** (`03_bed_to_bigwig/`) — The normalised bedgraph is converted to the binary indexed bigWig format using `bedGraphToBigWig`. These are the **recommended files for loading into IGV, UCSC Genome Browser, or for any cross-sample comparisons**.

3. **Downsampled bigWig** (`03_bed_to_bigwig_downsampled/`) — For purely visual comparison, BAMs are downsampled to the same read depth before conversion, eliminating sequencing-depth differences.

4. **IgG-subtracted bigWig** (`07_bigwig_minus_igg/`) — The IgG control bigWig is subtracted from each target sample bigWig using `bigwigCompare`, producing a background-corrected signal track.

### Peak Calling

#### SEACR (`04_called_peaks/seacr/`)
[SEACR](https://github.com/FredHutch/SEACR) (Sparse Enrichment Analysis for CUT&RUN) is specifically designed for low-background assays. It identifies regions with signal above either an empirical threshold (top X% of AUC; `--seacr_peak_threshold`) or above the IgG control bedgraph (recommended). SEACR output BED files contain:
- Columns 1–3: chromosome, start, end
- Column 4: Total signal (AUC) in the peak region
- Column 5: Max signal within the peak
- Column 6: Region of max signal

#### MACS2 (`04_called_peaks/macs2/`)
[MACS2](https://github.com/macs3-project/MACS) is a widely used peak caller that models the shift size and local background. Recommended for data with higher background noise or for transcription-factor narrow peaks. Key parameters exposed in this pipeline:
- `--macs2_pvalue` — p-value cutoff (default: 1e-5)
- `--macs2_gsize` — effective genome size
- `--macs2_narrow_peak` — use narrow-peak mode (for TFs; broad mode is default for histone marks)

### Consensus Peaks (`05_consensus_peaks/`)
Per-sample peaks for each experimental group are merged using `bedtools merge` to produce a set of consensus peak regions representing the union of peaks called in any replicate of the group. A replicate threshold (`--replicate_threshold`) can be applied to retain only peaks present in ≥ N replicates, increasing stringency.

### FRiP Score (`07_peak_qc/frip/`)
The **Fraction of Reads in Peaks** (FRiP) score is calculated as:

```
FRiP = fragments overlapping peaks / total fragments
```

Fragment BED files extracted from BAMs (`06_fragments_from_bams/`) are intersected with the called peaks using `bedtools intersect`. FRiP scores are summarised in MultiQC. **FRiP ≥ 0.3** is generally considered acceptable; high-quality experiments often reach ≥ 0.7.

### Merged Peaks Table (`08_merged_peaks_table/`)
All consensus peaks across all groups are merged into a single genomic coordinate set. Read counts from each sample are then quantified within this unified set using `bedtools coverage`, producing a matrix of peak × sample counts. This table is suitable for downstream differential binding analysis (e.g. DESeq2 or DiffBind).

### Homer Motif Analysis (`09_homer_motifs/`)
Enabled with `--run_homer_motifs true`. [Homer](http://homer.ucsd.edu/homer/motif/) performs hypergeometric-optimised motif enrichment on:
- **merged_peaks/** — the union of all consensus peaks
- **consensus_peaks/{group}_motifs/** — group-specific consensus peaks

Key output files per peak set:
| File | Description |
|------|-------------|
| `homerResults.html` | Interactive HTML with de-novo motif logos and statistics |
| `knownResults.txt` | Known motif enrichment table (name, consensus, p-value, % target peaks) |
| `homerMotifs.all.motifs` | Position weight matrices for all discovered motifs |
| `Known_Motifs_Summary.txt` | Cross-group comparison of top known motifs |
| `DeNovo_Motifs_Summary.txt` | Cross-group comparison of top de-novo motifs |

The motif search window size is set by `--homer_motif_size` (default: 200 bp). Use 50–200 for sharp histone marks (H3K4me3) or 500–1000 for broad marks (H3K27me3).

### Peak Feature Annotation (`10_peak_feature_annotation/`)
Enabled with `--run_homer_peak_annotation true`. Each peak is annotated with its nearest genomic feature (promoter, exon, intron, intergenic, etc.) using Homer `annotatePeaks.pl`. Outputs include:

**Raw Homer annotation** (`annotatepeaks_raw/`)
- `*.annotatePeaks.txt` — per-sample annotation tables from Homer, one row per peak, with columns for chromosome, start, end, peak name, score, strand, and genomic feature annotation (including distance to TSS and GC content via the `-CpG` flag).

**Aggregated summary** (`10_peak_feature_annotation/`)
Produced by the `SUMMARIZE_PEAK_ANNOTATIONS` step, which collects all per-sample Homer tables and generates:

| File | Description |
|------|-------------|
| `homer_peak_annotation.raw_annotation_summary.tsv` | Raw Homer annotation categories across all samples |
| `homer_peak_annotation.feature_summary.tsv` | Canonicalised feature categories (Promoter, Exon, Intron, Intergenic, etc.) |
| `homer_peak_annotation.sample_stats.tsv` | Per-sample peak counts and mean GC% |
| `homer_peak_annotation.feature_percent_table.tsv` | Percentage of peaks in each feature per sample |
| `homer_peak_annotation.gc_per_peak.tsv` | GC% for every individual peak |
| `homer_peak_annotation.gc_by_sample.tsv` | Sample-level GC% summary |
| `homer_peak_annotation.gc_by_sample.{png,pdf}` | GC distribution plot per sample |
| `homer_peak_annotation.functional_enrichment.tsv` | Gene Ontology Biological Process enrichment |
| `homer_peak_annotation.stacked_bar.{png,pdf}` | Stacked bar chart of feature composition across samples |

### ChIPseeker Enhanced Visualisations (`10_peak_feature_annotation/`)
Enabled with `--run_homer_peak_annotation true` **and** `--run_chipseeker true`. [ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html) is an R/Bioconductor package that reads the Homer annotation tables and produces a richer set of visualisations.

**Per-sample plots** (produced by `CHIPSEEKER_ANNOTATE`; one set per replicate):

| File | Description |
|------|-------------|
| `{sample}_chipseeker_annotation.tsv` | Full ChIPseeker annotation table for the sample |
| `{sample}_chipseeker_annotation_pie.{png,pdf}` | Pie chart of genomic feature distribution |
| `{sample}_chipseeker_annotation_bar.{png,pdf}` | Bar plot of genomic feature counts |
| `{sample}_chipseeker_annotation_upset.{png,pdf}` | Upset plot showing overlapping annotation categories |
| `{sample}_chipseeker_dist_to_tss.{png,pdf}` | Distribution of peak distances relative to the TSS |
| `{sample}_chipseeker_vennpie.{png,pdf}` | Venn-pie representation of annotation overlap |
| `{sample}_chipseeker_coverage.{png,pdf}` | Chromosome coverage plot |

**Cross-sample comparative plots** (produced by `CHIPSEEKER_COMPARE`; requires `--run_chipseeker_compare true`):

| File | Description |
|------|-------------|
| `chipseeker_comparison_annotation_bar.{png,pdf}` | Stacked annotation bar chart across all samples |
| `chipseeker_comparison_annotation_bar_slim_version.{png,pdf}` | Simplified annotation bar chart |
| `chipseeker_comparison_annotation_bar_condition.{png,pdf}` | Annotation bar chart grouped by experimental condition |
| `chipseeker_comparison_dist_to_tss.{png,pdf}` | TSS distance distribution comparison across all samples |
| `chipseeker_comparison_coverage.{png,pdf}` | Genome-wide coverage comparison |
| `chipseeker_comparison_coverage_condition_average.{png,pdf}` | Per-condition average coverage |
| `chipseeker_condition_{group}_pie.{png,pdf}` | Per-condition annotation pie charts |

**How ChIPseeker derives annotations:**
1. Homer's `annotatePeaks.pl` output is read and converted to GRanges objects.
2. `annotatePeak()` from ChIPseeker maps each peak to its nearest genomic feature using the TxDb annotation derived from the pipeline GTF.
3. Features are collapsed into canonical categories: Promoter (TSS ± 3 kb), Exon, Intron, Downstream, 5' UTR, 3' UTR, and Intergenic.
4. Comparative plots aggregate per-sample annotations to enable cross-group comparisons.

**Key parameters:**
- `--run_chipseeker` (default: false) — enable per-sample ChIPseeker annotation and plots
- `--run_chipseeker_compare` (default: true) — enable cross-sample comparative plots (only active when `run_chipseeker` is also true)

## Quality Thresholds (Recommended)

| Metric | Recommendation |
|--------|---------------|
| FRiP score | ≥ 0.3 (abundant marks); ≥ 0.1 (rare marks / TFs) |
| Peak count consistency across replicates | < 2× fold difference between biological replicates |
| Peak reproducibility | ≥ 50% of peaks reproduced in at least one other replicate |
