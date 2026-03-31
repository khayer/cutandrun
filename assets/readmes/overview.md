# CUT&RUN Pipeline Results — Overview

This directory contains the complete output of the **nf-core/cutandrun** pipeline (Weitzman Lab custom build).

## Result Subdirectories

| Directory | Contents |
|-----------|----------|
| [`00_genome/`](00_genome/README.md) | Reference genome FASTA, annotation, and Bowtie2 indices used for alignment |
| [`01_prealign/`](01_prealign/README.md) | Pre-alignment quality control: raw FastQC, adapter trimming (TrimGalore), and post-trim FastQC |
| [`02_alignment/`](02_alignment/README.md) | Bowtie2 alignments to target and spike-in genomes; duplicate marking; filtered BAM files |
| [`03_peak_calling/`](03_peak_calling/README.md) | Coverage tracks (bedgraph, bigWig), called peaks (SEACR/MACS2), consensus peaks, motif analysis, and peak annotation |
| [`04_reporting/`](04_reporting/README.md) | Quality-control visualisations: deepTools QC (PCA, fingerprint, correlation), heatmaps, upset plots, IGV session, pyGenomeTracks snapshots, and MultiQC summary |
| [`pipeline_info/`](pipeline_info/README.md) | Nextflow execution reports, software versions, and run parameters |

## Suggested Reading Order

1. **`04_reporting/multiqc/`** — Start here. The MultiQC report aggregates the most important QC metrics in a single interactive HTML page.
2. **`01_prealign/`** — Inspect raw read quality and trimming effectiveness.
3. **`02_alignment/`** — Evaluate alignment rates for target and spike-in genomes; review duplication metrics.
4. **`03_peak_calling/`** — Examine called peaks, consensus peaks, bigWig tracks, and motif results.
5. **`04_reporting/`** — Review deepTools QC, heatmaps, and genome-browser sessions for visual validation.

## Normalisation Note

Spike-in normalisation scales each sample's read coverage by the fraction of reads that aligned to the spike-in genome (E. coli or other exogenous reference), so that samples with different epitope abundances can be compared directly. If `--normalisation_mode_dual true` was used, a second round of within-experiment normalisation is also applied. The normalised bigWig tracks in `03_peak_calling/03_bed_to_bigwig/` should be used for all downstream comparisons.

## Key Parameters Used

| Parameter | Purpose |
|-----------|---------|
| `--normalisation_mode` | Spike-in, RPKM, or CPM normalisation |
| `--peakcaller` | `seacr` (default) or `macs2` |
| `--run_homer_motifs` | Enable de-novo and known motif discovery |
| `--run_peak_signal_profiler` | Enable PeakSignalProfiler module |
| `--run_pygenometracks_top10` | Enable locus snapshots of top consensus peaks |

For a full list of parameters and their meanings see [`pipeline_info/params.json`](pipeline_info/params.json) and the [pipeline documentation](https://nf-co.re/cutandrun/output).
