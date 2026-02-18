PeakSignalProfiler Nextflow wrapper

This workflow runs PeakSignalProfiler's `02_run_multisample.R` inside a Singularity image.

Usage example:

nextflow run workflows/peak_signal_profiler/main.nf \
  --samplesheet samplesheet_processed.csv \
  --annotation HSV17_from_gff.bed.bed.gz \
  --genome 17_No_repeats.fasta.fai \
  --psp_sif /path/to/peaksignalprofiler.sif \
  --psp_dir /path/to/peak-signal-profiler \
  -profile singularity -with-report -with-trace -resume

Notes:
- The process maps Nextflow CPUs to the R script via `--cores ${task.cpus}`.
- The R script is instructed to write into `psp_out/` inside the task workdir; that directory is published by Nextflow.
- Provide a valid Singularity image via `--psp_sif` (placeholder used by default).
- Optionally run multiple samplesheets by creating a channel and invoking the subworkflow `PEAK_SIGNAL_PROFILER` directly in a higher-level Nextflow script.

## Preflight package check ✅

This module performs a preflight inside the task/container to verify the following R/Bioconductor packages are available before running `02_run_multisample.R`:

- rtracklayer
- GenomicRanges
- IRanges
- S4Vectors
- ggplot2
- optparse

If any package is missing the process will fail early with a clear message instructing you to rebuild the Singularity image (SIF). Example error printed to stderr:

> ERROR: Required R packages missing inside container. Rebuild SIF to include: rtracklayer,GenomicRanges,IRanges,S4Vectors,ggplot2,optparse

How to verify packages inside a SIF:

```bash
singularity exec /path/to/psp-1.0.0.sif \
  Rscript -e "pkgs <- c('rtracklayer','GenomicRanges','IRanges','S4Vectors','ggplot2','optparse'); missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; if(length(missing)) { cat('MISSING:', paste(missing, collapse=','), '\n'); quit(status=2) } else { cat('OK\n') }"
```

How to add the packages (Singularity / Docker recipe snippet):

```
# inside %post of a Singularity definition or Dockerfile
Rscript -e "install.packages(c('ggplot2','optparse'), repos='https://cloud.r-project.org')"
Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('rtracklayer','GenomicRanges','IRanges','S4Vectors'))"
```

After rebuilding the SIF re-run the workflow with `--psp_sif /path/to/psp-1.0.0.sif`.
