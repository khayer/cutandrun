# DESeq2 Peak Differential Analysis Integration

This Nextflow pipeline now includes integrated DESeq2 differential abundance analysis for comparing peak enrichment between conditions.

## Quick Start

### 1. Prepare Your Contrasts File

Create a TSV file specifying which conditions to compare. Example: `deseq2_contrasts.txt`

```
treatment_name	control_name	treatment_col_indices	control_col_indices
RAD51_B02_26	RAD51_DMSO_26	17,18,19	22,23,24
PAA_TI_26	RAD51_DMSO_26	13,14	22,23,24
DRB_RI_26	RAD51_DMSO_26	8,9	22,23,24
```

**How to determine column indices:**
1. After peak calling and merging, review `results/03_peak_calling/08_merged_peaks_table/all_samples.tab`
2. The first 3 columns are chromosome, start, end
3. Starting from column 4, each subsequent column is a sample's read counts
4. Count from column 4 onward (1-based indexing where column 1 = chr, 2 = start, 3 = end, 4 = first sample, etc.)
5. Example: If "RAD51_B02_26_R1" is the 14th column in the table, use column index 14

### 2. Run the Pipeline

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --genome HSV17 \
  --run_deseq2_peak_analysis \
  --deseq2_contrasts deseq2_contrasts.txt \
  -profile {singularity,conda,docker}
```

Or add to your nextflow.config:

```nextflow
params {
    run_deseq2_peak_analysis = true
    deseq2_contrasts = "deseq2_contrasts.txt"
}
```

Then run simply: `nextflow run main.nf -profile singularity`

### 3. Access Results

Results are in: `results/03_peak_calling/09_deseq2_analysis/`

For each contrast (e.g., `RAD51_B02_26_vs_RAD51_DMSO_26/`):
- `RAD51_B02_26_vs_RAD51_DMSO_26_results.csv` - Full results for all peaks
- `RAD51_B02_26_vs_RAD51_DMSO_26_significant.csv` - Peaks with padj < 0.05
- `RAD51_B02_26_vs_RAD51_DMSO_26_significant_log2fc.bedGraph` - BedGraph of significant peaks (for IGV)
- `RAD51_B02_26_vs_RAD51_DMSO_26_volcano_plot.pdf` - Volcano plot with top 5 genes labeled
- `RAD51_B02_26_vs_RAD51_DMSO_26_ma_plot.pdf` - MA plot

## Contrasts File Format

**Required Format:** Tab-separated values (TSV)

**Columns:**
1. `treatment_name` - Name of treatment condition (e.g., "RAD51_B02_26")
2. `control_name` - Name of control condition (e.g., "RAD51_DMSO_26")
3. `treatment_col_indices` - Comma-separated column indices for treatment replicates
4. `control_col_indices` - Comma-separated column indices for control replicates

**Notes:**
- Column indices are 1-based and refer to columns in `all_samples.tab` starting from column 4 (first sample)
- Use commas (no spaces) to separate multiple column indices: `17,18,19`
- Lines starting with `#` are treated as comments
- Empty lines are skipped
- Each row defines one contrast

**Example contrasts file:**

```
# HSV-17 virus genome peak differential analysis
# Format: treatment_name control_name treatment_cols control_cols

RAD51_B02_26	RAD51_DMSO_26	17,18,19	22,23,24
PAA_TI_26	RAD51_DMSO_26	13,14	22,23,24
DRB_RI_26	RAD51_DMSO_26	8,9	22,23,24
```

## How to Determine Column Indices

1. **Locate the counts table:**
   ```bash
   results/03_peak_calling/08_merged_peaks_table/all_samples.tab
   ```

2. **View the header (first line):**
   ```bash
   head -1 results/03_peak_calling/08_merged_peaks_table/all_samples.tab
   ```

3. **Map your samples to columns:**
   - Columns 1-3: `#chr`, `start`, `end`
   - Column 4 onward: Sample counts (one per sample, in sorted order)
   - Count from column 4 as index 4, column 5 as index 5, etc.

4. **Find your sample's indices:**
   For example, if `all_samples.tab` has:
   ```
   #chr    start   end     ADNP_B02_25_R1  ADNP_DMSO_25_R1 ... RAD51_B02_26_R1 RAD51_B02_26_R2 RAD51_B02_26_R3 RAD51_DMSO_26_R1 ...
   HSV17   100     200     ...
   ```
   Then RAD51_B02_26 replicates are at columns 17, 18, 19 and RAD51_DMSO_26 controls are at columns 22, 23, 24.

## Output Files Explained

### `*_results.csv`
Complete DESeq2 results table for all peaks analyzed:
- `peak_id`: Peak coordinates (chr:start:end)
- `chrom`: Chromosome
- `start`, `end`: Peak boundaries
- `nearest_gene`: Nearest annotated gene
- `nearest_gene_distance`: Distance to nearest gene (bp)
- `display_label`: Gene name or peak ID for visualization
- `log2FoldChange`: log2(treatment/control) abundance
- `lfcSE`: Standard error of log2FC
- `stat`: Test statistic
- `pvalue`: Raw p-value
- `padj`: Adjusted p-value (Benjamini-Hochberg FDR)

### `*_significant.csv`
Filtered results containing only peaks with **padj < 0.05** (statistically significant)

### `*_significant_log2fc.bedGraph`
4-column BED format for IGV visualization:
```
HSV17   1000    1200    1.234
HSV17   2000    2100   -0.567
```
Where the last column is log2FoldChange for each significant peak. Load into IGV to visualize peak changes across the genome.

### Volcano Plot PDF
Scatter plot of log2FoldChange vs. -log10(p-value):
- **X-axis**: log2 Fold Change (left=depleted, right=enriched in treatment)
- **Y-axis**: -log10(p-value) (higher = more significant)
- **Red points**: Significant peaks (padj < 0.05)
- **Black points**: Non-significant peaks
- **Labeled points**: Top 5 most significant peaks in each direction
- **Blue dashed lines**: Significance thresholds (padj=0.05, log2FC=±1)

### MA Plot PDF
Relationship between peak abundance and magnitude of change:
- **X-axis**: log2(mean counts) - average abundance
- **Y-axis**: log2 Fold Change
- **Red points**: Significant peaks
- Shows if differential abundance is driven by highly or lowly abundant peaks

## Analysis Details

### Peak Selection
- Union of all peaks called in either treatment or control group
- Ensures inclusion of both enhanced and depleted regions
- Avoids bias from only analyzing constitutively present peaks

### Differential Testing
- **Method**: DESeq2 (Love et al., 2014)
- **Statistics**: 
  - Size factor normalization (built-in to DESeq2)
  - Negative binomial model with dispersion estimation
  - Wald test for condition contrast
  - Multiple testing correction: Benjamini-Hochberg FDR
- **Significance threshold**: padj < 0.05 (standard for genomics)

### Gene Annotation
Each peak is annotated with:
- Nearest gene (from GTF file)
- Distance to that gene (base pairs)
- Gene name used in volcano plot labels

## Troubleshooting

### Error: "No narrowPeak files found"
- Check that peak calling completed successfully
- Verify sample names match between `all_samples.tab` and narrowPeak files
- Expected path: `results/03_peak_calling/04_called_peaks/macs2/{sample_name}.macs2_peaks.narrowPeak`

### Error: "Invalid column indices"
- Verify indices are integers, not sample names
- Check that indices are comma-separated without spaces: `17,18,19` ✓ not `17, 18, 19`
- Ensure indices match actual columns in `all_samples.tab`

### No significant peaks found
- May be expected if conditions are very similar
- Check raw results CSV to see log2FC distribution
- Verify correct samples were specified in contrasts file

### Different results than standalone R script
- May occur if different genome annotation (GTF) is used
- Check that all_samples.tab is identical
- Verify sample column indices are correct

## Integration with Other Analyses

DESeq2 peak analysis runs **after**:
- Peak calling (MACS2/SEACR)
- Peak consensus/merging
- BAM file quality control

And **alongside**:
- Promoter GC differential binding analysis (if enabled)
- ChIPseeker annotation
- Other QC steps

Results are independent and can be run together or separately.

## Next Steps After DESeq2 Analysis

1. **Filter results** by your desired log2FC threshold (e.g., |log2FC| > 1)
2. **Load bedGraph** files into IGV for genome-wide visualization
3. **Validate findings** by examining raw CUT&RUN signal tracks
4. **Cross-reference genes** with functional annotations or literature
5. **Perform pathway analysis** on significant genes (using nearest gene annotation)
6. **Integrate with RNA-seq** data if available (compare peak changes with transcript changes)

## References

- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
- MACS2: Zhang Y, et al. (2008). Model-based Analysis of ChIP-Seq (MACS). *Genome Biology*, 9(9), R137.
