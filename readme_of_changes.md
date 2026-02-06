# Custom Modifications to nf-core/cutandrun Pipeline

This document tracks custom modifications made to the nf-core/cutandrun pipeline for the Weitzman Lab.

## Table of Contents
- [Dual Normalization Feature](#dual-normalization-feature)
- [BigWig Subtraction](#bigwig-subtraction)
- [Merged Peaks Table](#merged-peaks-table)
- [Read Count Annotation](#read-count-annotation)

---

## Dual Normalization Feature

**Date Added:** February 6, 2026  
**Purpose:** Enable normalization by both spike-in genome (E. coli) and target genome (e.g., HSV) abundance

### Background
Standard spike-in normalization corrects for technical variation in ChIP efficiency between samples using:
```
scale_factor = normalisation_c / spike_in_aligned_reads
```

However, this doesn't account for biological variation in target genome abundance (e.g., different viral loads, genome copy numbers). The dual normalization feature addresses both:

```
scale_factor = (normalisation_c / spike_in_reads) × (mean_target_reads / sample_target_reads)
```

This provides:
- **Technical normalization** (spike-in): Corrects for ChIP efficiency differences
- **Biological normalization** (target genome): Corrects for viral load or genome abundance differences

### Usage

Enable dual normalization by adding the parameter:

```bash
nextflow run main.nf \
  --normalisation_mode Spikein \
  --normalisation_mode_dual true \
  --normalisation_c 10000 \
  ... other parameters ...
```

**Requirements:**
- Must be used with `--normalisation_mode Spikein`
- Requires spike-in genome alignment (default: E. coli K12-MG1655)
- Target genome alignment statistics must be available

### Implementation Details

**Files Modified:**

1. **`nextflow.config`**
   - Added parameter: `normalisation_mode_dual = false` (default: disabled)

2. **`subworkflows/local/prepare_peakcalling.nf`**
   - Added inputs: `target_metadata` and `mean_target_reads`
   - Modified spike-in normalization to:
     - Load both spike-in and target genome alignment metadata
     - Calculate mean target reads across all samples
     - Apply dual factor: `spikein_factor × target_abundance_factor`
   - Added logging to show individual sample scale factors

3. **`workflows/cutandrun.nf`**
   - Added channel to calculate mean target genome reads
   - Passes `ch_metadata_bt2_target` and `ch_mean_target_reads` to PREPARE_PEAKCALLING
   - Calculation only performed when `normalisation_mode_dual = true`

4. **`nextflow_schema.json`**
   - Added `normalisation_mode_dual` parameter definition with help text

### Output Changes

When enabled, scale factors in bigWig files and bedGraph files will reflect both normalizations. The pipeline will log each sample's:
- Spike-in normalization factor
- Target genome abundance factor
- Final combined scale factor

### Example Use Cases

1. **Viral CUT&RUN experiments:** Normalize for both ChIP efficiency and viral genome copy number across infection time points
2. **Copy number variation studies:** Account for both technical and biological variation when genome abundance varies
3. **Multi-condition experiments:** Compare samples with different target genome abundance (e.g., drug treatments affecting viral replication)

### Notes

- If target genome alignment data is missing or zero, falls back to spike-in normalization only
- Mean target reads calculated only from samples with valid alignment counts (> 0)
- Dual normalization factor is multiplicative, not additive
- Compatible with all downstream peak calling and visualization steps

---

## BigWig Subtraction

**Date Added:** Prior to February 6, 2026  
**Purpose:** Generate background-subtracted bigWig files (ChIP - IgG control)

### Usage

```bash
nextflow run main.nf \
  --run_bigwig_subtract true \
  --use_control true \
  ... other parameters ...
```

### Files Generated

Located in: `results/03_peak_calling/07_bigwig_minus_igg/`

For each target sample:
- `{sample}.log2ratio.bigWig` - Log2(ChIP/Control)
- `{sample}.subtract.bigWig` - ChIP - Control (linear subtraction)

### Implementation

- Uses `deeptools bigwigCompare` with both `--operation log2` and `--operation subtract`
- Control matching: Targets with replicate N use control replicate N
- Channel matching: Uses `.combine()` + `.filter()` for many-to-one mapping (multiple targets can use same control)
- Integrated into IGV session XML automatically

---

## Merged Peaks Table

**Date Added:** Prior to February 6, 2026  
**Purpose:** Create Homer mergePeaks-style presence/absence matrix for peaks across all samples

### Usage

Automatically runs when peak calling is enabled.

### Output

Located in: `results/03_peak_calling/08_merged_peaks_table/`

- `merged_peaks_table.txt` - Tab-delimited table with peak presence (1/0) per sample
- `merged_peaks.bed` - BED file of merged peak regions

### Table Format

```
PeakID  Chr    Start    End      Length  Sample1  Sample2  Sample3  Total
Peak_1  chr1   1000     1500     500     1        0        1        2
Peak_2  chr1   2000     2300     300     1        1        1        3
```

### Implementation

- Module: `modules/local/python/merge_peaks_table.nf`
- Merges overlapping peaks using 50% overlap threshold
- Creates binary presence/absence matrix
- Counts total samples with each peak

---

## Read Count Annotation

**Date Added:** Prior to February 6, 2026  
**Purpose:** Annotate merged peaks with read counts from all BAM files using deeptools multiBamSummary

### Usage

Automatically runs after merged peaks table generation.

### Output

Located in: `results/03_peak_calling/08_merged_peaks_table/`

- `all_samples.npz` - Binary numpy matrix (for deeptools)
- `all_samples.tab` - Tab-delimited read counts per peak per sample

### Table Format

```
#chr   start   end     Sample1.bam     Sample2.bam     Sample3.bam
chr1   1000    1500    245             12              189
chr1   2000    2300    567             523             601
```

### Implementation

- Module: `modules/nf-core/deeptools/multibamsummary_bed/main.nf`
- Uses `multiBamSummary BED-file` mode
- Counts reads overlapping each merged peak region
- Includes all BAM files (targets and controls)
- Sample order is sorted alphabetically for consistency

### Use Cases

- Differential binding analysis
- Quantitative comparison of peak intensities
- Statistical testing between conditions
- Input for downstream clustering/heatmaps

---

## Version Information

**Pipeline Base:** nf-core/cutandrun v3.2.2  
**Nextflow Version:** >=23.04.0  
**Custom Modifications Tested With:**
- Singularity containers
- SLURM executor
- HSV-1 viral genome (136 kb)
- E. coli K12-MG1655 spike-in

---

## Validation Status

All custom modifications have been tested with:
- Sample data: 11 target samples + 3 IgG controls
- Genome: HSV-1 strain 17 (136 kb)
- Peak caller: MACS2
- Execution: Successfully completed end-to-end pipeline runs

---

## Contact

**Lab:** Weitzman Lab  
**Modified By:** GitHub Copilot (AI Assistant) - Claude Sonnet 4.5
**User:** hayerk  
**Institution:** CHOP

For questions or issues with these modifications, please refer to the commit history or contact the lab directly.
