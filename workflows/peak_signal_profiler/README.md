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
