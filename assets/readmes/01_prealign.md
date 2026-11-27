# 01_prealign — Pre-alignment Quality Control

This directory contains quality-control data for raw and trimmed sequencing reads, generated **before** any genome alignment is performed.

## Subdirectories

| Path | Contents |
|------|----------|
| `merged_fastq/` | Concatenated FASTQ files (only when `--save_merged_fastq` is set and multiple runs were provided per sample) |
| `pretrim_fastqc/` | FastQC reports on the **raw**, untrimmed reads |
| `trimgalore/` | Adapter-trimmed FASTQ files (only when `--save_trimmed` is set) and TrimGalore log files |
| `trimgalore/fastqc/` | FastQC reports on the **trimmed** reads |

## How These Files Are Derived

### FASTQ Merging (`merged_fastq/`)
When multiple sequencing runs are provided for the same sample (e.g. for read-depth top-up), the pipeline concatenates all FASTQ files for that sample before any further processing. Enable saving with `--save_merged_fastq`.

### Raw FastQC (`pretrim_fastqc/`)
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is run on the raw (or merged) input reads. Each HTML report contains:
- **Sequence counts** — total reads per sample
- **Per-base quality scores** — should remain ≥ Q30 for modern Illumina data
- **Per-sequence GC content** — should match the expected GC% of the target organism
- **Overrepresented sequences / adapter content** — adapter contamination is expected here and will be removed by TrimGalore

> **CUT&RUN note:** Some FastQC modules (e.g. per-base sequence content) will flag as "FAIL" for CUT&RUN data. This is normal and is caused by Tn5 preferential binding or the 10-bp sawtooth periodicity characteristic of the assay. These warnings do **not** indicate a failed experiment.

### Adapter Trimming (`trimgalore/`)
[TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) wraps [Cutadapt](https://cutadapt.readthedocs.io/) and FastQC to remove adapter sequences and low-quality bases. Adapters are detected automatically. The `*_trimming_report.txt` log files record:
- Which adapter was detected and trimmed
- Number and percentage of reads that were quality-trimmed
- Number of reads that were discarded (too short after trimming)

### Post-trim FastQC (`trimgalore/fastqc/`)
FastQC is re-run on the trimmed reads. Compare these reports with the pre-trim reports to confirm that:
- Adapter content has been removed
- Overrepresented sequences are reduced (any remaining may be biologically relevant)
- Read quality remains acceptable

## Quality Thresholds (Recommended)

| Metric | Recommendation |
|--------|---------------|
| Raw read count | ≥ 6 M per sample (for abundant marks like H3K27me3; more for low-abundance targets) |
| Post-trim Q30 | ≥ Q30 across the majority of read length |
| Adapter content (post-trim) | < 1% at any position |
