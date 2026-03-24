# ChIPseeker Integration Summary - nf-core/cutandrun Workflow

## 1. Overall Workflow Structure

### Main Entry Point
- **File**: `main.nf` → delegates to `workflows/cutandrun.nf` via `CUTANDRUN()` workflow
- **Pattern**: DSL2 (Nextflow v2+) with modular includes

### Major Workflow Stages (in `workflows/cutandrun.nf`)

1. **Genome Preparation** (`PREPARE_GENOME`)
   - Outputs: FASTA, GTF, bowtie2 indices, spike-in genomes
   - Conditional: `params.run_genome_prep`

2. **Input Check** (`INPUT_CHECK`)
   - Validates sample sheet
   - Outputs: standardized fastq channels with metadata
   - Channel structure: `[[id, group, replicate, single_end, is_control], [reads]]`

3. **Read QC & Trimming** (`FASTQC_TRIMGALORE`)
   - FastQC and TrimGalore processing
   - Conditional: `params.run_trim_galore_fastqc`

4. **Alignment** (`ALIGN_BOWTIE2`)
   - Aligns to both target and spike-in genomes (if `normalisation_mode == "Spikein"`)
   - Outputs: BAM files, alignment logs
   - Conditional: `params.run_alignment`

5. **Read Filtering & Deduplication**
   - Filter reads: `FILTER_READS` subworkflow
   - Mark/remove duplicates: `MARK_DUPLICATES_PICARD`
   - Optional linear dedup: `DEDUPLICATE_LINEAR`

6. **Peak Calling Preparation** (`PREPARE_PEAKCALLING`)
   - Converts BAM → bedgraph/bigwig with normalization
   - Supports dual normalization (spike-in + target abundance)

7. **Peak Calling** (MACS2 or SEACR)
   - Multiple peak callers supported
   - Generates primary peaks per replicate

8. **Consensus Peaks** (`CONSENSUS_PEAKS`)
   - Merges replicate peaks by group
   - Conditional: `params.run_consensus_peaks`

9. **Peak Annotation Branch** ← **Primary ChIPseeker Integration Point**
   - See Section 2 below

10. **Optional Analyses**
    - Homer motif discovery: `run_homer_motifs`
    - Peak signal profiler: `run_peak_signal_profiler`
    - Heatmaps (DeepTools)
    - IGV session generation

---

## 2. Peak Annotation Workflow Details

### Current Implementation (HOMER-based)

**Location in workflow**: Lines 796-825 in `workflows/cutandrun.nf`

```
Input: ch_peaks_primary (replicate-level peak BEDs)
  ↓
HOMER_ANNOTATEPEAKS (if run_homer_peak_annotation == true)
  - For each replicate: annotatePeaks.pl with GTF, CpG calculation
  - Output: {sample}.annotatePeaks.txt (tab-separated table)
  ↓
SUMMARIZE_PEAK_ANNOTATIONS (if run_homer_peak_annotation == true)
  - Aggregates all *.annotatePeaks.txt files
  - Generates summary tables and stacked bar plots
  - Output: homer_peak_annotation.*
```

### Channel Flow for Peak Annotation

**Input Channel** (`ch_peaks_primary`):
```groovy
[[id: sample_id, group: group_name, replicate: rep_num], bed_file]
```

**HOMER_ANNOTATEPEAKS Process**:
- Input: tuple of (meta, bed), fasta file, gtf file
- Output: tuple of (meta, annotated_table)
- Uses direct `-CpG` flag for GC% calculation

**SUMMARIZE_PEAK_ANNOTATIONS Module**:
- Input: Collection of all *.annotatePeaks.txt files
- Outputs (11 files):
  - `homer_peak_annotation.raw_annotation_summary.tsv` → raw HOMER categories
  - `homer_peak_annotation.feature_summary.tsv` → canonicalized features
  - `homer_peak_annotation.sample_stats.tsv` → per-sample peak counts & GC%
  - `homer_peak_annotation.feature_percent_table.tsv` → normalized percentages
  - `homer_peak_annotation.gc_per_peak.tsv` → individual peak GC%
  - `homer_peak_annotation.gc_by_sample.tsv` → sample-level GC summary
  - `homer_peak_annotation.gc_by_sample.{png,pdf}` → GC distribution plots
  - `homer_peak_annotation.functional_enrichment.tsv` → GO BP enrichment
  - `homer_peak_annotation.stacked_bar.{png,pdf}` → feature composition plots

### Output Directories
- Raw HOMER outputs: `results/03_peak_calling/10_peak_feature_annotation/annotatepeaks_raw/`
- Summary reports: `results/03_peak_calling/10_peak_feature_annotation/`

---

## 3. Module File Structure Patterns

### Directory Organization
```
modules/local/
├── linux/              # AWK, CUT, other shell tools
├── python/             # Python analysis modules
│   ├── summarize_peak_annotations.nf
│   ├── merge_peaks_table.nf
│   ├── frag_len_hist.nf
│   └── ...
├── homer/              # HOMER-specific tools
│   ├── annotatepeaks/main.nf
│   └── findmotifsgenome/main.nf
└── samtools_custom_view.nf
```

### Module Template Pattern (from `homer/annotatepeaks/main.nf`)

```groovy
process MODULE_NAME {
    tag "$meta.id"                    // Identifies each process execution
    label 'process_medium'            // Resource allocation (from configs)

    conda "bioconda::tool=version"    // Conda environment
    container "...singularity..."     // Container URI

    input:
    tuple val(meta), path(bed)        // Standard: tuple(metadata, files)
    path fasta
    path gtf

    output:
    tuple val(meta), path("${prefix}.annotatePeaks.txt"), emit: annot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when    // Conditional execution

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    command here...
    """

    stub:                             // For testing/dry-run
    // Minimal stub output
}
```

### Naming Conventions

| Element | Pattern | Example |
|---------|---------|---------|
| Process name | ALL_CAPS_WITH_UNDERSCORES | `HOMER_ANNOTATEPEAKS` |
| Module file | lowercase_with_underscores.nf | `summarize_peak_annotations.nf` |
| Emit label | lowercase_with_underscores | `emit: raw_summary` |
| Variable prefix | camelCase | `ch_peaks_primary`, `ch_software_versions` |
| Output file | snake_case_with_identifiers | `{prefix}.annotatePeaks.txt` |

### Key Metadata Structure

All modules expect `meta` tuple:
```groovy
meta = [
    id: "sample_name",           // Unique sample identifier
    group: "group_name",         // Experimental group
    replicate: 1,                // Replicate number
    single_end: false,           // Read type
    is_control: false            // Control flag
]
```

---

## 4. Current `summarize_peak_annotations.py` Details

### Location
- Binary: `bin/summarize_peak_annotations.py`
- Module wrapper: `modules/local/python/summarize_peak_annotations.nf`

### Input Processing
```python
# Reads all *.annotatePeaks.txt files generated by HOMER
# Expected column names (case-insensitive):
# - "Annotation" or "annotation" (feature assignment)
# - "GC%", "gc", or similar (GC content %)
```

### Feature Canonicalization
HOMER annotations → 8 standardized feature classes:
- Promoter (includes TSS annotations)
- Exon
- Intron
- 5' UTR / 3' UTR
- Intergenic
- Non-coding
- TTS (Transcription Termination Site)
- Unknown

**Function**: `canonical_feature(annotation: str) → str`

### Key Processing Functions

1. **`parse_homer_table(path)`**
   - Returns: (raw_df, feature_df, gc_per_peak_df, gene_list, stats_dict)
   - Extracts gene names from HOMER's "Nearest Gene" columns
   - Handles GC% as fractions or percentages
   - Counts peaks per feature category

2. **`extract_gene_list(df)`**
   - Parses HOMER's gene annotation columns (can have multiple genes/peak)
   - Cleans tokens separated by commas, semicolons, slashes

3. **`run_functional_enrichment(gene_sets)`**
   - Uses gseapy for Gene Ontology Brain Biological Process enrichment
   - Returns top terms per sample

4. **`main()`**
   - Aggregates all samples into unified tables
   - Generates stacked bar plot with N and GC labels per replicate
   - Outputs 11+ analysis files

### Output Specifications

| Output File | Format | Purpose |
|-------------|--------|---------|
| `raw_annotation_summary.tsv` | TSV | HOMER's native annotation categories |
| `feature_summary.tsv` | TSV | Canonicalized features summary |
| `feature_percent_table.tsv` | TSV | Feature percentages (for plotting) |
| `sample_stats.tsv` | TSV | N peaks, mean GC% per replicate |
| `gc_per_peak.tsv` | TSV | Individual peak GC% values |
| `gc_by_sample.tsv` | TSV | Sample-level GC statistics |
| `functional_enrichment.tsv` | TSV | GO BP enrichment terms |
| `stacked_bar.{png,pdf}` | Image | Main visualization |
| `gc_by_sample.{png,pdf}` | Image | GC distribution plots |

---

## 5. Parameter Naming Conventions

### Workflow Control Parameters

**Pattern**: `run_*`, `skip_*`, `only_*` prefixes

| Parameter | Type | Default | Stage | Purpose |
|-----------|------|---------|-------|---------|
| `run_genome_prep` | boolean | true | Genome setup | Skip/run genome file prep |
| `run_input_check` | boolean | true | Input validation | Skip/run samplesheet check |
| `run_alignment` | boolean | true | Alignment | Skip/run read alignment |
| `run_mark_dups` | boolean | true | QC | Skip/run deduplication |
| `run_homer_peak_annotation` | boolean | false | Peak annotation | Enable HOMER annotatePeaks |
| `run_homer_motifs` | boolean | false | Motif discovery | Enable Homer motif finding |
| `run_peak_signal_profiler` | boolean | false | Analysis | Skip/run peak profiler |
| `skip_fastqc` | boolean | false | QC | Skip FastQC |

### Tool-Specific Parameters

**HOMER Peak Annotation**:
- `homer_peak_annotation_tss_dist` (default: "1000,2000,3000")
  - Passed directly to `annotatePeaks.pl -d` flag
  - Controls TSS/promoter distance windows (comma-separated)

**HOMER Motif Finding**:
- `run_homer_motifs` (boolean, default: false)
- `homer_motif_size` (int, default: 200)
  - Peak region size for motif discovery (bp)
  - Common values: 50-200 (sharp marks), 500-1000 (broad marks), 200 (promoters)

**Normalization**:
- `normalisation_mode` (default: "Spikein")
  - Options: "Spikein", "RC", "Spike-in+Target" (dual normalization)
- `igg_scale_factor` (default: 0.5)
  - Scaling factor for IgG control samples

**Peak Calling**:
- `peakcaller` (default: "macs2")
  - Options: "macs2", "seacr", or both (comma-separated)
- `use_control` (boolean, default: true)
  - Use control samples in peak calling

**Consensus Peaks**:
- `consensus_peak_mode` (default: "group")
  - Which peaks to merge (per-group or all-merge)
- `replicate_threshold` (default: 1.0)
  - Fraction of replicates that must have peak for inclusion

### Global Output Control

- `outdir` (default: "./results")
  - Output directory
- `skip_reporting` (boolean, default: false)
  - Skip all reporting modules
- `skip_multiqc` (boolean, default: false)
  - Skip MultiQC report generation
- `publish_dir_mode` (default: "copy")
  - How to handle output files (copy, link, move, etc.)

---

## 6. Integration Points for ChIPseeker

### Recommended Locations

1. **Parallel Implementation** (recommended):
   - Create `modules/local/chipseeker/annotatepeaks/main.nf`
   - Add parameter: `run_chipseeker_peak_annotation` (boolean)
   - Add to peak annotation conditional block (line ~796)
   - Process flow identical to HOMER but using R/ChIPseeker

2. **Alternative: Replace HOMER**
   - Modify existing `HOMER_ANNOTATEPEAKS` module
   - Add parameter: `peak_annotation_tool` (values: "homer", "chipseeker")
   - Conditional logic to select tool

3. **Key Requirements**:
   - Accept input: tuple(meta, bed), fasta, gtf
   - Output: tuple(meta, annotation_table)
   - Format output as TSV with standardized columns
   - Compatible with `SUMMARIZE_PEAK_ANNOTATIONS` or custom summarizer

### Module Discovery
- Local modules: checked in `modules/local/` first
- Included in workflows via: `include { MODULE_NAME } from "path/to/main.nf"`
- All modules require explicit includes (no auto-discovery)

