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
| `10_peak_feature_annotation/` | Genomic feature annotation of peaks (only when `--run_homer_peak_annotation true`) |

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
- Raw annotation tables (`annotatepeaks_raw/`)
- Summary bar charts and tables of feature distribution across groups

## Quality Thresholds (Recommended)

| Metric | Recommendation |
|--------|---------------|
| FRiP score | ≥ 0.3 (abundant marks); ≥ 0.1 (rare marks / TFs) |
| Peak count consistency across replicates | < 2× fold difference between biological replicates |
| Peak reproducibility | ≥ 50% of peaks reproduced in at least one other replicate |
