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
    Rscript -e "pkgs <- c('rtracklayer','GenomicRanges','IRanges','S4Vectors','ggplot2','optparse','cowplot'); missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; if(length(missing)) { cat('MISSING:', paste(missing, collapse=','), '\n'); quit(status=2) } else { cat('OK\n') }"


CI / common failure (building slim SIF)

- Symptom: GitHub Actions step that builds a slim SIF fails with the message:

  "glib-2.0 headers are required to build common."

- Cause: Apptainer/Singularity requires system development headers (glib, libgpgme, libseccomp, etc.) on the runner when building from source.

- Quick fixes:
  1. In your GitHub Actions job (Ubuntu runner) **install the build deps before compiling Apptainer**. Example step to add to the workflow job that installs Apptainer:

```yaml
- name: Install Apptainer build deps
  run: |
    sudo apt-get update
    sudo apt-get install -y build-essential automake autoconf libtool m4 pkg-config \
      libseccomp-dev libgpgme-dev libssl-dev uuid-dev squashfs-tools libarchive-dev \
      libglib2.0-dev libcap-dev
```

  2. Or run the helper script included in this repo (below) from your workflow.

- Rationale: `libglib2.0-dev` supplies the missing glib-2.0 headers; the other packages are commonly required to build Apptainer successfully.


Helper: install Apptainer build deps (use in CI)

  ./dev/docker/peaksignalprofiler/install_apptainer_deps.sh


