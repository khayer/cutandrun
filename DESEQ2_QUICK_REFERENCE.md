# DESeq2 Peak Analysis - Quick Reference

## Command-Line Quick Start

```bash
# 1. Create contrasts.txt file (see example below)
# 2. Run pipeline
nextflow run main.nf \
  --input samplesheet.csv \
  --genome HSV17 \
  --run_deseq2_peak_analysis \
  --deseq2_contrasts contrasts.txt \
  -profile singularity

# 3. View results
ls -la results/03_peak_calling/09_deseq2_analysis/
```

## Contrasts File Format (TSV)

```
treatment_name        control_name          treatment_cols      control_cols
RAD51_B02_26          RAD51_DMSO_26         17,18,19            22,23,24
PAA_TI_26             RAD51_DMSO_26         13,14               22,23,24
DRB_RI_26             RAD51_DMSO_26         8,9                 22,23,24
```

## How to Find Column Indices

```bash
# View header of all_samples.tab
head -1 results/03_peak_calling/08_merged_peaks_table/all_samples.tab | tr '\t' '\n' | nl

# Count columns manually and note sample columns
# Columns 1-3 are: chr, start, end
# Column 4 onward are sample counts
# Use these column numbers in contrasts file
```

## Pipeline Configuration

Add to `nextflow.config`:
```nextflow
params {
    run_deseq2_peak_analysis = true
    deseq2_contrasts = "contrasts.txt"
}
```

Then run: `nextflow run main.nf -profile singularity`

## Output Interpretation

### Results CSV Columns
- `log2FoldChange`: Enrichment in treatment (positive = more in treatment)
- `pvalue`: Raw p-value from Wald test
- `padj`: Adjusted p-value (use this for significance, threshold: 0.05)
- `nearest_gene`: Annotated gene nearest to peak
- `nearest_gene_distance`: Distance in base pairs

### Key Outputs
| File | Use |
|------|-----|
| `*_results.csv` | All peaks (for filtering) |
| `*_significant.csv` | Peaks with padj < 0.05 |
| `*_significant_log2fc.bedGraph` | IGV visualization |
| `*_volcano_plot.pdf` | Overview of results |

## Common Tasks

### Filter for enriched peaks (log2FC > 1)
```R
res <- read.csv("results.csv")
enriched <- res[res$log2FoldChange > 1 & res$padj < 0.05, ]
```

### Load bedGraph in IGV
1. Open IGV
2. File → Load from File
3. Select `*_significant_log2fc.bedGraph`
4. View alongside BAM files

### Check column mapping
```bash
awk 'NR==1 {for(i=1;i<=NF;i++) print i": "$i}' \
  results/03_peak_calling/08_merged_peaks_table/all_samples.tab | head -30
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| No narrowPeak files found | Check sample names match between samplesheet and peak files |
| Invalid column indices | Use 1-based integers, comma-separated, no spaces: `17,18,19` |
| DESeq2 estimation failed | May occur with very few peaks; check that consensus peak calling succeeded |
| No significant peaks | Expected if conditions very similar; check log2FC distribution in results CSV |

## File Locations

```
Pipeline Config:    nextflow.config
Module:             modules/local/r/deseq2_peak_analysis/main.nf
R Script:           bin/deseq2_peak_analysis.R
Documentation:      docs/DESEQ2_PEAK_ANALYSIS.md
Example Contrasts:  assets/deseq2_contrasts.txt
```

## Key Parameters

```nextflow
run_deseq2_peak_analysis = false        # Set to true to enable
deseq2_contrasts = null                 # Path to contrasts TSV file
outdir = "./results"                    # Output directory (from main config)
```

## Interpretation Guide

**log2 Fold Change**
- `> 0` = Higher in treatment
- `< 0` = Higher in control
- `|log2FC| > 1` = At least 2-fold difference

**Adjusted P-value (padj)**
- `< 0.05` = Significant (FDR-corrected)
- `> 0.05` = Not significant at standard threshold

**Volcano Plot**
- **Red points** = Significant (padj < 0.05)
- **X-axis crossing** = log2FC = 1 (2-fold threshold)
- **Y-axis height** = Significance level

## Default Peak Selection

- **Union of Peaks**: All peaks called in treatment OR control group
- **Rationale**: Captures both enhanced and depleted regions
- **Not**: Intersection (only peaks in both groups)

## Performance Notes

- Processing time: 2-5 minutes per contrast (depends on peak count)
- Memory: ~4 GB per analysis (can be adjusted with `label 'process_medium'`)
- Storage: ~5 MB per contrast output

## References

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

---

**Full documentation**: See `docs/DESEQ2_PEAK_ANALYSIS.md`  
**Integration details**: See `DESEQ2_INTEGRATION.md`
