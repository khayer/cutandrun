#!/usr/bin/env Rscript
# Install ChIPseeker and dependencies if not already installed

required_packages <- c(
    "ChIPseeker",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "ggplot2",
    "GenomicRanges",
    "GenomicFeatures"
)

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install Bioconductor packages if needed
if (!require("BiocManager")) install.packages("BiocManager")

for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
        if (pkg %in% c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", 
                       "GenomicRanges", "GenomicFeatures")) {
            BiocManager::install(pkg, ask = FALSE)
        } else {
            install.packages(pkg)
        }
    }
}

cat("All ChIPseeker dependencies installed successfully!\n")
