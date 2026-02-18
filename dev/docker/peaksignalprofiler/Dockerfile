# Dockerfile for PeakSignalProfiler runtime image
# Builds an R image with the Bioconductor/CRAN packages required by
# PeakSignalProfiler (rtracklayer, GenomicRanges, IRanges, S4Vectors, ggplot2, optparse)

FROM rocker/r-ver:4.3.1

LABEL maintainer="Katharina Hayer <khayer@weitzmanlab.org>"
LABEL org.opencontainers.image.source="https://github.com/hayerk/peak-signal-profiler"

ENV DEBIAN_FRONTEND=noninteractive

# system deps commonly required by Bioconductor packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    && rm -rf /var/lib/apt/lists/*

# install minimal CRAN + Bioconductor packages the module checks for
RUN Rscript -e "install.packages(c('ggplot2','optparse'), repos='https://cloud.r-project.org')"
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('rtracklayer','GenomicRanges','IRanges','S4Vectors'))"

# runtime working dir (pipeline will mount scripts via --psp_dir)
WORKDIR /opt/psp

# image metadata
ARG BUILD_DATE
ARG VCS_REF
LABEL org.opencontainers.image.created=$BUILD_DATE
LABEL org.opencontainers.image.revision=$VCS_REF

CMD ["Rscript", "--version"]
