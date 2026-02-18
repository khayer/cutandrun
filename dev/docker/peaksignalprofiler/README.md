Build & publish (CI) — PeakSignalProfiler image

Quick notes
- CI builds a Docker image and pushes to GitHub Container Registry (GHCR).
- To use Singularity (.sif) build locally from the published Docker image (example below).

Local SIF build (example):

  # build Docker locally (optional)
  docker build -t peaksignalprofiler:local -f dev/docker/peaksignalprofiler/Dockerfile .

  # build SIF from Docker Hub / GHCR image (requires Singularity/Apptainer locally)
  singularity build psp-1.0.0.sif docker://ghcr.io/<OWNER>/peaksignalprofiler:latest

Verify R packages inside SIF:

  singularity exec psp-1.0.0.sif \
    Rscript -e "pkgs <- c('rtracklayer','GenomicRanges','IRanges','S4Vectors','ggplot2','optparse'); missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; if(length(missing)) { cat('MISSING:', paste(missing, collapse=','), '\n'); quit(status=2) } else { cat('OK\n') }"
