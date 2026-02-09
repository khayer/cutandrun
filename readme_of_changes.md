# Custom Modifications to nf-core/cutandrun Pipeline

**Author:** Katharina Hayer (khayer)  
**Co-created with:** GitHub Copilot (Claude Sonnet 4.5)  

This document tracks custom modifications made to the nf-core/cutandrun pipeline for the Weitzman Lab.

## How I use it:

    ./nextflow run main.nf   --input samplesheet_just2026.csv   --outdir results_dual_norm   --normalisation_mode Spikein   --normalisation_mode_dual true   --normalisation_c 10000   -profile singularity   --fasta HSV17_genome_files/17_No_repeats.fasta   --bowtie2 /home/hayerk/data/index/weitzman_human_hsv/HSV17_genome/   --gene_bed HSV17_from_gff.bed   --gtf HSV17_genome_files/17_No_repeats.gff   --macs_gsize 136000   -work-dir work_dual_norm   --peakcaller macs2 -resume --run_homer_motifs true --skip_multiqc --bigwigcompare_binsize 5

    ./nextflow run main.nf   --input samplesheet.csv   --outdir results_human   --normalisation_mode Spikein   --normalisation_mode_dual true   --normalisation_c 10000   -profile singularity    -work-dir work_human   --peakcaller macs2 --run_homer_motifs true --skip_multiqc --bigwigcompare_binsize 5 --genome hg38


## Table of Contents
- [Homer Motif Analysis](#homer-motif-analysis)
- [Dual Normalization Feature](#dual-normalization-feature)
- [BigWig Subtraction](#bigwig-subtraction)
- [Merged Peaks Table](#merged-peaks-table)
- [Read Count Annotation](#read-count-annotation)

---

## Homer Motif Analysis

**Date Added:** February 6, 2026  
**Purpose:** Perform de novo and known motif discovery in peak regions using Homer

### Background

[Homer](http://homer.ucsd.edu/homer/motif/) (Hypergeometric Optimization of Motif EnRichment) is a comprehensive suite for motif discovery and next-generation sequencing analysis. This integration adds motif finding capabilities to identify:
- **De novo motifs**: Novel transcription factor binding motifs enriched in peaks
- **Known motifs**: Enrichment of previously characterized motifs from databases
- **Motif locations**: Precise positions of motifs within peak regions

Motif analysis is performed on:
1. **Merged peaks** - Combined peak set across all samples
2. **Consensus peaks** - Group-specific consensus peaks (per experimental group)

### Usage

Enable Homer motif analysis:

```bash
nextflow run main.nf \
  --run_homer_motifs true \
  --homer_motif_size 200 \
  ... other parameters ...
```

**Parameters:**
- `--run_homer_motifs` (default: false) - Enable/disable Homer motif finding
- `--homer_motif_size` (default: 200) - Size of region for motif finding (bp)
  - Use 200 for promoters
  - Use `given` to use exact peak size
  - Use 50-200 for sharp histone marks (H3K4me3)
  - Use 500-1000 for broad marks (H3K27me3)

### Output

Located in: `results/03_peak_calling/09_homer_motifs/`

**Merged peaks:**
- `merged_peaks/merged_peaks_motifs/` - Complete Homer output directory
- `merged_peaks/merged_peaks_motifs/homerResults.html` - Main results page
- `merged_peaks/merged_peaks_motifs/knownResults.txt` - Known motif enrichment table
- `merged_peaks/merged_peaks_motifs/homerMotifs.all.motifs` - All discovered de novo motifs

**Consensus peaks (per group):**
- `consensus_peaks/{group}_motifs/` - One directory per experimental group
- `consensus_peaks/{group}_motifs/homerResults.html` - Group-specific results
- `consensus_peaks/{group}_motifs/knownResults.txt` - Group-specific known motifs

**Summary Reports:**
- `Known_Motifs_Summary.txt` - Comprehensive comparison of known motifs across all groups
- `DeNovo_Motifs_Summary.txt` - Comprehensive comparison of de novo motifs across all groups

### Key Output Files

| File | Description |
|------|-------------|
| `homerResults.html` | Interactive HTML with motif logos, statistics, and target sequences |
| `knownResults.txt` | Table of known motif enrichment (p-values, % targets, % background) |
| `motifN.motif` | Individual de novo motif files (position weight matrices) |
| `homerMotifs.all.motifs` | Combined file of all discovered motifs |
| `motifN.similar.motifs.html` | Similar known motifs for each de novo motif |
| **`Known_Motifs_Summary.txt`** | **Cross-group comparison of top 10 known motifs** |
| **`DeNovo_Motifs_Summary.txt`** | **Cross-group comparison of top 10 de novo motifs** |

### Summary Reports Format

The automatically generated summary reports provide:

**Per-Group Analysis:**
- Top 10 motifs ranked by p-value
- Motif name, consensus sequence, p-value, and % of peaks containing the motif
- Separate sections for merged peaks and each consensus peak group

**Cross-Group Comparison:**
- Shared motifs between groups
- Unique motifs per group
- Side-by-side comparison of enrichment statistics
- Merged peaks vs. consensus peaks comparison

**Example Summary Output:**
```
====================================================================================================
Peak Set: merged_peaks
====================================================================================================

  Rank   Motif                                    Consensus             P-value          % Target  
  ------ ---------------------------------------- -------------------- --------------- ----------
  1      RAD51                                    GCTGGGCG             1e-89            45.3      
  2      E2F1                                     TTTCGCGC             1e-42            28.7      
  3      TP53                                     RRRCWWGYYY           1e-23            18.2      

====================================================================================================
MOTIF COMPARISON ACROSS GROUPS
====================================================================================================

DRB_RI_26 vs PAA_TI_26:
  Shared motifs: 8
    - RAD51
    - E2F1
    - TP53
  Unique to DRB_RI_26: 2
  Unique to PAA_TI_26: 2
```

### Interpreting Results

**Known Motifs Table (`knownResults.txt`):**
- **Motif Name**: Transcription factor or motif ID
- **Consensus**: Best matching sequence
- **P-value**: Statistical significance of enrichment
- **Log P-value**: -log10(p-value) for visualization
- **% of Target**: Percentage of peaks containing motif
- **% of Background**: Percentage in background sequences

**De novo Motifs:**
- Ranked by enrichment p-value
- Include sequence logos and position weight matrices
- Can be compared against known motif databases

### Implementation Details

**Files Modified:**

1. **`modules/local/homer/findmotifsgenome/main.nf`**
   - New module wrapping `findMotifsGenome.pl`
   - Inputs: BED file, genome FASTA, motif size
   - Outputs: Complete Homer directory, HTML results, known motifs table

2. **`workflows/cutandrun.nf`**
   - Added `HOMER_FINDMOTIFSGENOME_MERGED` for merged peaks analysis
   - Added `HOMER_FINDMOTIFSGENOME_CONSENSUS` for consensus peaks analysis
   - Integrated after peak merging/consensus steps

3. **`conf/modules.config`**
   - Added publishDir: `03_peak_calling/09_homer_motifs/`
   - Configured CPU usage: `-p 4` (4 cores for motif finding)
   - Separate outputs for merged and consensus peaks
   - Summary reports published to top-level motifs directory

4. **`nextflow.config`**
   - Added parameters: `run_homer_motifs` and `homer_motif_size`

5. **`modules/local/python/summarize_homer_motifs.nf`**
   - New module that parses all Homer results
   - Generates Known_Motifs_Summary.txt and DeNovo_Motifs_Summary.txt
   - Compares motifs across groups automatically
   - Runs after all Homer analyses complete

### Use Cases

1. **Transcription factor discovery:** Identify enriched TF binding sites in ChIP-seq peaks
2. **Co-factor analysis:** Find motifs co-occurring with primary target
3. **Condition comparison:** Compare motifs between different experimental groups
4. **Motif evolution:** Track motif changes across time points or treatments
5. **Validation:** Confirm expected TF binding in immunoprecipitation experiments

### Notes

- Homer automatically selects background regions matched for GC content
- Motif finding is computationally intensive; runs with 4 CPUs by default
- Results are best interpreted in context of known biology/literature
- Multiple testing correction is applied automatically
- De novo motifs are compared against JASPAR, TRANSFAC, and other databases

### Example Output Interpretation

For a RAD51 ChIP-seq experiment:
```
Known Motifs:
  1. RAD51(Homeo)/K562-RAD51-ChIP-Seq  P-value: 1e-50  % Targets: 45.2%
  2. E2F1(E2F)/Hela-E2F1-ChIP-Seq      P-value: 1e-23  % Targets: 28.7%
```

This suggests:
- Primary RAD51 binding motif is highly enriched (expected)
- E2F1 co-factor binding is also enriched (interesting biological insight)

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
- Spike-in: E. coli K12-MG1655
- Execution: Successfully completed end-to-end pipeline runs
- Homer: v4.11 (motif analysis)

---

## Contact

**Lab:** Weitzman Lab  
**Modified By:** GitHub Copilot (AI Assistant) - Claude Sonnet 4.5
**User:** hayerk  
**Institution:** CHOP

For questions or issues with these modifications, please refer to the commit history or contact the lab directly.
