# DESeq2 Peak Differential Analysis - Nextflow Integration Summary

## What Was Added

The DESeq2 differential peak analysis has been integrated as a permanent module in the Nextflow CUT&RUN pipeline.

### New Files Created

1. **Nextflow Module**: `modules/local/r/deseq2_peak_analysis/main.nf`
   - Defines the DESEQ2_PEAK_ANALYSIS process
   - Handles input/output channels
   - Manages conda/container environments
   - Publishes results to dedicated output directory

2. **R Analysis Script**: `bin/deseq2_peak_analysis.R`
   - Standalone script called by the Nextflow module
   - Accepts command-line arguments for flexibility
   - Performs complete differential abundance analysis
   - Generates publication-quality plots
   - Produces bedGraph files for genome visualization

3. **Documentation**: `docs/DESEQ2_PEAK_ANALYSIS.md`
   - Comprehensive user guide
   - Format specifications for contrasts file
   - Example contrasts and workflows
   - Troubleshooting guide
   - Output file explanations

4. **Example Contrasts File**: `assets/deseq2_contrasts.txt`
   - Ready-to-use example with HSV-17 samples
   - Well-commented template for user customization

### Modified Files

1. **nextflow.config**
   - Added `run_deseq2_peak_analysis` flag (default: false)
   - Added `deseq2_contrasts` parameter (path to contrasts TSV file)

2. **workflows/cutandrun.nf**
   - Imported DESEQ2_PEAK_ANALYSIS module
   - Added workflow logic to parse contrasts file
   - Integrated analysis into main pipeline execution
   - Added software version tracking

## How to Use

### Basic Usage

1. **Create contrasts file** (TSV format):
   ```
   treatment_name	control_name	treatment_col_indices	control_col_indices
   RAD51_B02_26	RAD51_DMSO_26	17,18,19	22,23,24
   ```

2. **Run pipeline**:
   ```bash
   nextflow run main.nf \
     --input samplesheet.csv \
     --genome HSV17 \
     --run_deseq2_peak_analysis \
     --deseq2_contrasts deseq2_contrasts.txt \
     -profile singularity
   ```

3. **Access results**:
   ```
   results/03_peak_calling/09_deseq2_analysis/{treatment}_vs_{control}/
   ```

### Output Files per Contrast

- `*_results.csv` - Full DESeq2 results (all peaks)
- `*_significant.csv` - Significant peaks only (padj < 0.05)
- `*_significant_log2fc.bedGraph` - BedGraph for IGV visualization
- `*_volcano_plot.pdf` - Volcano plot with gene labels
- `*_ma_plot.pdf` - MA plot

## Key Features

✅ **Flexible Contrasts**: Specify any treatment vs. control comparison  
✅ **Multiple Comparisons**: Run multiple contrasts in single pipeline execution  
✅ **Gene Annotation**: Automatic nearest-gene annotation for each peak  
✅ **Publication-Ready Plots**: Volcano and MA plots with significance highlighting  
✅ **IGV Compatible**: bedGraph output for direct genome visualization  
✅ **Error Handling**: Comprehensive error messages for common issues  
✅ **Reproducible**: Full containerization and conda environments specified  

## Technical Details

### Differential Testing Method

- **Statistical Framework**: DESeq2 (Love et al., 2014)
- **Normalization**: DESeq2's built-in size factor estimation
- **Model**: Negative binomial generalized linear model
- **Test**: Wald test for condition contrast
- **Multiple Testing Correction**: Benjamini-Hochberg FDR

### Peak Selection Strategy

- **Union of Consensus Peaks**: All peaks called in either treatment OR control group
- **Rationale**: Captures both gain and loss of signal
- **Avoids Bias**: Doesn't restrict to peaks only called in one group

### Column Index Determination

The contrasts file uses 1-based column indices into the `all_samples.tab` file:
- Columns 1-3: chr, start, end (always present)
- Columns 4+: Sample read counts (sorted alphabetically by sample name)

Example:
```
#chr    start   end     ADNP_B02_25_R1  ...  RAD51_B02_26_R1  R2  R3  RAD51_DMSO_26_R1  R2  R3
                         4               ...  17               18  19  22                 23  24
```

## Integration Points

### Inputs
- **Peak narrowPeak files**: `results/03_peak_calling/04_called_peaks/macs2/*.narrowPeak`
- **Merged counts table**: `results/03_peak_calling/08_merged_peaks_table/all_samples.tab`
- **GTF annotation**: From PREPARE_GENOME module or user-specified
- **Contrasts file**: User-provided TSV with contrast definitions

### Dependencies
- Runs **after**: Peak calling and peak consensus/merging (automatic)
- Runs **alongside**: Other analyses (ChIPseeker, promoter GC, etc.)
- Can be toggled independently via `--run_deseq2_peak_analysis`

### Output Directory
```
results/03_peak_calling/09_deseq2_analysis/
├── {treatment}_vs_{control}/
│   ├── {treatment}_vs_{control}_results.csv
│   ├── {treatment}_vs_{control}_significant.csv
│   ├── {treatment}_vs_{control}_significant_log2fc.bedGraph
│   ├── {treatment}_vs_{control}_volcano_plot.pdf
│   └── {treatment}_vs_{control}_ma_plot.pdf
```

## Default Parameters

```nextflow
params {
    run_deseq2_peak_analysis = false    // Disabled by default
    deseq2_contrasts = null              // User must provide
}
```

## Reusable Design Notes

### From Standalone Script → Pipeline Integration

The original standalone R script (`deseq2_RAD51_comparison.R`) was refactored for pipeline integration:

1. **Parameterization**: Hard-coded contrasts → command-line arguments
2. **Flexibility**: Three fixed contrasts → unlimited user-specified contrasts
3. **Reproducibility**: Manual execution → automated within pipeline
4. **Containerization**: System R dependencies → conda/container specification
5. **Error Handling**: Silent failures → informative error messages
6. **Result Management**: Manual output organization → automatic publishDir

### Backward Compatibility

The standalone script remains available and unchanged:
```bash
module load R/4.4.0
Rscript deseq2_RAD51_comparison.R
```

Both the standalone and pipeline-integrated versions can be used independently.

## Quality Assurance

✅ Module follows nf-core best practices  
✅ All dependencies specified in conda/container  
✅ Error handling for missing files and invalid parameters  
✅ Software versions tracked in outputs  
✅ Comprehensive logging and user feedback  
✅ Example contrasts file provided  
✅ Full documentation included  

## Future Enhancements (Optional)

Potential improvements for future versions:
- Support for interaction terms (more complex designs)
- Volcano plot filtering thresholds as parameters
- Custom color schemes for plots
- Export to common formats (GCT, GFF3)
- Integration with pathway analysis tools
- Summary report generation

## Support

For detailed usage information, see: `docs/DESEQ2_PEAK_ANALYSIS.md`

For troubleshooting, check the documentation's "Troubleshooting" section or review the Nextflow execution logs.
