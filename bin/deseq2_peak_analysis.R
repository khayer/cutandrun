#!/usr/bin/env Rscript

# DESeq2 Peak Differential Abundance Analysis
# For use in Nextflow pipeline
# Requires: DESeq2, GenomicRanges, tidyverse, rtracklayer, optparse

suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(GenomicRanges)
    library(rtracklayer)
    library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--counts-table"), type = "character", help = "Path to merged peaks count table"),
  make_option(c("--peak-dir"), type = "character", help = "Path to directory with narrowPeak files"),
  make_option(c("--gtf-file"), type = "character", help = "Path to GTF annotation file"),
  make_option(c("--samplesheet"), type = "character", help = "Path to samplesheet CSV file"),
  make_option(c("--treatment-group"), type = "character", help = "Treatment group name from samplesheet"),
  make_option(c("--control-group"), type = "character", help = "Control group name from samplesheet"),
  make_option(c("--output-prefix"), type = "character", help = "Prefix for output files")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# Read samplesheet to find samples in each group
samplesheet <- read.csv(args$samplesheet, stringsAsFactors = FALSE)
colnames(samplesheet) <- tolower(colnames(samplesheet))

treatment_rows <- samplesheet[samplesheet$group == args$`treatment-group`, ]
control_rows <- samplesheet[samplesheet$group == args$`control-group`, ]

if (nrow(treatment_rows) == 0) {
  stop(paste("No samples found for treatment group:", args$`treatment-group`))
}
if (nrow(control_rows) == 0) {
  stop(paste("No samples found for control group:", args$`control-group`))
}

# Create sample IDs: groupname_R{replicate}
treatment_samples <- paste0(treatment_rows$group, "_R", treatment_rows$replicate)
control_samples <- paste0(control_rows$group, "_R", control_rows$replicate)

cat("Treatment samples found:", paste(treatment_samples, collapse=", "), "\n")
cat("Control samples found:", paste(control_samples, collapse=", "), "\n")

# Read all_samples.tab header to find column indices
counts_header <- readLines(args$`counts-table`, n = 1)
header_cols <- strsplit(counts_header, "\t")[[1]]

# Find column indices for our samples (columns 1-3 are chr, start, end)
treatment_indices <- which(header_cols %in% treatment_samples)
control_indices <- which(header_cols %in% control_samples)

if (length(treatment_indices) == 0) {
  cat("Available columns:", paste(header_cols[4:min(10, length(header_cols))], collapse=", "), "...\n")
  stop(paste("Could not find columns for treatment samples:", paste(treatment_samples, collapse=", ")))
}
if (length(control_indices) == 0) {
  cat("Available columns:", paste(header_cols[4:min(10, length(header_cols))], collapse=", "), "...\n")
  stop(paste("Could not find columns for control samples:", paste(control_samples, collapse=", ")))
}

cat("Treatment column indices:", paste(treatment_indices, collapse=", "), "\n")
cat("Control column indices:", paste(control_indices, collapse=", "), "\n")

# Read narrowPeak files
read_narrowpeak <- function(sample_name, peak_dir) {
  file_path <- file.path(peak_dir, paste0(sample_name, ".macs2_peaks.narrowPeak"))
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
  GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2] + 1, df[, 3]))
}

# Read peak union for a set of samples
read_peak_union <- function(sample_names, peak_dir) {
  # Extract base sample name (remove _R{replicate})
  base_names <- gsub("_R[0-9]+$", "", sample_names)
  peak_list <- lapply(unique(base_names), function(x) read_narrowpeak(x, peak_dir))
  peak_list <- peak_list[!sapply(peak_list, is.null)]
  if (length(peak_list) == 0) {
    stop("No narrowPeak files were loaded for: ", paste(unique(base_names), collapse = ", "))
  }
  IRanges::reduce(Reduce(c, peak_list))
}

# Get gene label from GTF features
get_gene_label <- function(gene_features) {
  gene_label <- mcols(gene_features)$gene
  if (is.null(gene_label)) {
    gene_label <- mcols(gene_features)$Name
  }
  gene_label <- as.character(gene_label)
  fallback_name <- as.character(mcols(gene_features)$Name)
  fallback_locus <- as.character(mcols(gene_features)$locus_tag)
  fallback_id <- as.character(mcols(gene_features)$ID)
  gene_label[is.na(gene_label) | gene_label == ""] <- fallback_name[is.na(gene_label) | gene_label == ""]
  gene_label[is.na(gene_label) | gene_label == ""] <- fallback_locus[is.na(gene_label) | gene_label == ""]
  gene_label[is.na(gene_label) | gene_label == ""] <- fallback_id[is.na(gene_label) | gene_label == ""]
  gene_label
}

# Annotate peaks with nearest genes
annotate_peaks_with_genes <- function(res_df, chrom_filt, start_filt, end_filt, peak_ids, gtf_file) {
  peak_ranges <- GRanges(
    seqnames = sub("^HSV17$", "HSV", as.character(chrom_filt)),
    ranges = IRanges(start_filt + 1, end_filt),
    peak_id = peak_ids
  )

  gtf <- rtracklayer::import(gtf_file)
  gene_features <- gtf[mcols(gtf)$type == "gene"]
  gene_label <- get_gene_label(gene_features)

  nearest_gene_idx <- GenomicRanges::nearest(peak_ranges, gene_features)
  nearest_gene <- gene_label[nearest_gene_idx]
  nearest_distance <- GenomicRanges::distance(peak_ranges, gene_features[nearest_gene_idx])

  res_df$chrom <- chrom_filt
  res_df$start <- start_filt
  res_df$end <- end_filt
  res_df$nearest_gene <- nearest_gene
  res_df$nearest_gene_distance <- nearest_distance
  res_df$display_label <- ifelse(is.na(res_df$nearest_gene) | res_df$nearest_gene == "", res_df$peak_id, res_df$nearest_gene)

  res_df[, c(
    "peak_id", "chrom", "start", "end", "nearest_gene", "nearest_gene_distance",
    "display_label", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )]
}

# Write bedGraph for significant peaks
write_significant_bedgraph <- function(res_df, prefix) {
  sig_df <- res_df[!is.na(res_df$pvalue) & res_df$pvalue < 0.05, ]
  bedgraph_file <- paste0(prefix, "_significant_log2fc.bedGraph")
  if (nrow(sig_df) == 0) {
    file.create(bedgraph_file)
    return(invisible(bedgraph_file))
  }

  bedgraph_df <- data.frame(
    chrom = sig_df$chrom,
    start = pmax(sig_df$start - 1, 0),
    end = sig_df$end,
    log2FoldChange = sig_df$log2FoldChange
  )
  write.table(
    bedgraph_df,
    file = bedgraph_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  invisible(bedgraph_file)
}

# Main analysis
cat("\n=== Running DESeq2 Peak Analysis ===\n")
cat("Treatment:", args$`treatment-group`, "\n")
cat("Control:", args$`control-group`, "\n")

# Read peak unions
treatment_union <- read_peak_union(treatment_samples, args$`peak-dir`)
control_union <- read_peak_union(control_samples, args$`peak-dir`)

cat(args$`treatment-group`, "union peaks:", length(treatment_union), "\n")
cat(args$`control-group`, "union peaks:", length(control_union), "\n")

final_peaks <- IRanges::reduce(c(treatment_union, control_union))
cat("Peaks in either group (union):", length(final_peaks), "\n")

# Load count table
counts_df <- read.delim(args$`counts-table`, skip = 1)
chrom <- counts_df[, 1]
start <- counts_df[, 2]
end <- counts_df[, 3]
merged_gr <- GRanges(seqnames = as.character(chrom), ranges = IRanges(start + 1, end))

# Find overlapping peaks
overlaps <- findOverlaps(merged_gr, final_peaks)
keep_indices <- unique(queryHits(overlaps))
cat("Peaks from merged table matching consensus:", length(keep_indices), "\n")

counts_df_filtered <- counts_df[keep_indices, ]
chrom_filt <- counts_df_filtered[, 1]
start_filt <- counts_df_filtered[, 2]
end_filt <- counts_df_filtered[, 3]

# Extract counts for treatment and control samples
all_samples <- c(treatment_samples, control_samples)
all_indices <- c(treatment_indices, control_indices)

counts_subset <- as.data.frame(apply(counts_df_filtered[, all_indices], 2, as.integer))
colnames(counts_subset) <- all_samples
peak_ids <- paste(chrom_filt, start_filt, end_filt, sep = ":")
rownames(counts_subset) <- peak_ids
cat("Final peaks in DESeq2 analysis:", nrow(counts_subset), "\n")

# Create metadata
sample_metadata <- data.frame(
  sample = colnames(counts_subset),
  condition = c(
    rep(args$`treatment-group`, length(treatment_samples)),
    rep(args$`control-group`, length(control_samples))
  ),
  row.names = colnames(counts_subset)
)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_subset),
  colData = sample_metadata,
  design = ~ condition
)
dds$condition <- relevel(dds$condition, ref = args$`control-group`)

cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", args$`treatment-group`, args$`control-group`))

# Process results
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
res_df <- annotate_peaks_with_genes(res_df, chrom_filt, start_filt, end_filt, peak_ids, args$`gtf-file`)

# Write outputs
prefix <- args$`output-prefix`
write.csv(res_df, file = paste0(prefix, "_results.csv"), row.names = FALSE)

res_sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
write.csv(res_sig, file = paste0(prefix, "_significant.csv"), row.names = FALSE)

write_significant_bedgraph(res_df, prefix)

# Generate plots
pdf(paste0(prefix, "_ma_plot.pdf"), width = 8, height = 6)
plotMA(res, ylim = c(-3, 3))
dev.off()

pdf(paste0(prefix, "_volcano_plot.pdf"), width = 8, height = 6)
volcano_y <- -log10(res_df$pvalue + 1e-300)
plot(
  res_df$log2FoldChange,
  volcano_y,
  main = paste0("Volcano Plot: ", args$`treatment-group`, " vs ", args$`control-group`),
  xlab = "log2 Fold Change",
  ylab = "-log10(p-value)",
  pch = 19,
  cex = 0.7,
  col = ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "red", "black")
)
abline(v = c(-1, 1), h = -log10(0.05), lty = 2, col = "blue")

positive_df <- res_df[res_df$log2FoldChange > 0 & !is.na(res_df$pvalue), ]
positive_df <- positive_df[order(positive_df$pvalue), ]
positive_df <- head(positive_df, min(5, nrow(positive_df)))

negative_df <- res_df[res_df$log2FoldChange < 0 & !is.na(res_df$pvalue), ]
negative_df <- negative_df[order(negative_df$pvalue), ]
negative_df <- head(negative_df, min(5, nrow(negative_df)))

label_df <- unique(rbind(positive_df, negative_df))
if (nrow(label_df) > 0) {
  for (i in seq_len(nrow(label_df))) {
    point_x <- label_df$log2FoldChange[i]
    point_y <- -log10(label_df$pvalue[i] + 1e-300)
    point_label <- paste0(label_df$display_label[i], "\n", label_df$peak_id[i])
    text(
      point_x,
      point_y,
      labels = point_label,
      pos = ifelse(point_x >= 0, 4, 2),
      cex = 0.65,
      offset = 0.4
    )
  }
}

legend("topright", legend = "padj < 0.05", col = "red", pch = 19, bty = "n")
dev.off()

# Print summary
cat("\n=== DESeq2 Summary ===\n")
cat("Total peaks analyzed:", nrow(res_df), "\n")
cat("Significant peaks (padj < 0.05):", nrow(res_sig), "\n")
cat("Upregulated in", args$`treatment-group`, "(log2FC > 0):", sum(res_df$log2FoldChange > 0 & res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Downregulated in", args$`treatment-group`, "(log2FC < 0):", sum(res_df$log2FoldChange < 0 & res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Results saved with prefix:", prefix, "\n\n")
