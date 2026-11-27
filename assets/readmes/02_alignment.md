# 02_alignment — Alignment to Target and Spike-in Genomes

This directory contains BAM files and associated quality metrics produced by aligning trimmed reads to the target and spike-in reference genomes.

## Subdirectories

| Path | Contents |
|------|----------|
| `bowtie2/target/` | Initial sorted BAM files aligned to the target genome |
| `bowtie2/target/log/` | Bowtie2 alignment logs |
| `bowtie2/target/unmapped/` | Unmapped reads (only when `--save_unaligned` is set) |
| `bowtie2/target/markdup/` | BAM files after PCR duplicate **marking** (Picard MarkDuplicates) — these are the primary BAM files used downstream |
| `bowtie2/target/markdup/picard_metrics/` | Picard duplication metrics files |
| `bowtie2/target/dedup/` | BAM files after PCR duplicate **removal** (only when `--dedup_target_reads true`) |
| `bowtie2/target/linear_dedup/` | BAM files after Linear Amplification duplicate removal (only when applicable) |
| `bowtie2/spikein/` | BAM files aligned to the spike-in genome (only when `--save_spikein_aligned` is set) |
| `bowtie2/spikein/log/` | Spike-in Bowtie2 alignment logs |

## How These Files Are Derived

### 1. Bowtie2 Alignment
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) maps trimmed reads to both the target genome and the spike-in genome independently. Key settings:
- `--no-mixed` and `--no-discordant` flags enforce concordant paired-end alignment
- `--very-sensitive` mode is used for maximum sensitivity

Alignment rates for the **target** genome should be high (> 80% for a well-prepared library). Alignment rates for the **spike-in** genome are expected to be much lower (typically < 5% for target samples, 2–5% for IgG controls).

### 2. SAMtools Quality Filtering
BAM files are filtered using [SAMtools](http://www.htslib.org/) to:
- Remove reads below the minimum MAPQ threshold (`--minimum_alignment_q_score`, default: 20)
- Optionally remove reads aligned to mitochondrial chromosomes (`--remove_mitochondrial true`)
- Retain only properly paired reads

### 3. Picard MarkDuplicates (`markdup/`)
[Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) identifies reads that are likely PCR duplicates (same 5′ start position for both mates). Duplicates are **marked** but not removed by default, allowing downstream tools to optionally exclude them. IgG control samples are **additionally de-duplicated** (duplicates removed) as a separate step because the IgG is used as a control and lower background improves peak calling accuracy. Target sample de-duplication is off by default but can be enabled with `--dedup_target_reads true`.

The `*.MarkDuplicates.metrics.txt` files report:
- `ESTIMATED_LIBRARY_SIZE` — estimated number of unique molecules in the library
- `PERCENT_DUPLICATION` — fraction of reads identified as duplicates

> **Interpretation:** Moderate duplication levels (< 30%) are normal for CUT&RUN. Very high duplication (> 50%) in target samples may indicate over-amplification and can affect peak calling.

### 4. Linear Amplification Deduplication (`linear_dedup/`)
In assays using linear amplification, reads with the same read-1 5′ start site but different 3′ ends are expected biological duplicates. A custom deduplication step removes these by retaining only the highest-MAPQ read at each 5′ position.

### 5. Library Complexity (Preseq)
[Preseq](http://smithlabresearch.org/software/preseq/) models the relationship between sequencing depth and unique read yield to estimate library complexity. Results are reported in MultiQC (`04_reporting/multiqc/`). A rapidly saturating curve indicates the library was over-sequenced.

## Key Files

| File pattern | Description |
|---|---|
| `*.markdup.bam` | Primary BAM — used for all downstream peak calling and QC |
| `*.markdup.bam.bai` | BAI index for the primary BAM |
| `*.MarkDuplicates.metrics.txt` | Picard duplication statistics |
| `bowtie2/target/log/*.log` | Bowtie2 alignment rate summary |

## Quality Thresholds (Recommended)

| Metric | Recommendation |
|--------|---------------|
| Target alignment rate | ≥ 80% uniquely aligned |
| Spike-in alignment rate (IgG) | 2–10% (confirms spike-in was added) |
| % Duplication | < 30% for target samples |
| Aligned reads (post-filter) | ≥ 5 M for abundant marks |
