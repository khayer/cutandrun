# DESeq2 analysis for multiple group pairs against RAD51_DMSO_26

library(DESeq2)
library(tidyverse)
library(GenomicRanges)

peak_dir <- "results_HSV17_final/03_peak_calling/04_called_peaks/macs2"
counts_file <- "results_HSV17_final/03_peak_calling/08_merged_peaks_table/all_samples.tab"
gtf_file <- "HSV.gtf"
out_root <- "results_HSV17_final/deseq2_group_comparisons"

comparison_defs <- list(
  list(
    name = "RAD51_B02_vs_DMSO",
    treatment = "RAD51_B02_26",
    control = "RAD51_DMSO_26",
    treatment_samples = c("RAD51_B02_26_R1", "RAD51_B02_26_R2", "RAD51_B02_26_R3"),
    control_samples = c("RAD51_DMSO_26_R1", "RAD51_DMSO_26_R2", "RAD51_DMSO_26_R3"),
    treatment_indices = c(17, 18, 19),
    control_indices = c(22, 23, 24)
  ),
  list(
    name = "PAA_TI_26_vs_RAD51_DMSO_26",
    treatment = "PAA_TI_26",
    control = "RAD51_DMSO_26",
    treatment_samples = c("PAA_TI_26_R1", "PAA_TI_26_R2"),
    control_samples = c("RAD51_DMSO_26_R1", "RAD51_DMSO_26_R2", "RAD51_DMSO_26_R3"),
    treatment_indices = c(13, 14),
    control_indices = c(22, 23, 24)
  ),
  list(
    name = "DRB_RI_26_vs_RAD51_DMSO_26",
    treatment = "DRB_RI_26",
    control = "RAD51_DMSO_26",
    treatment_samples = c("DRB_RI_26_R1", "DRB_RI_26_R2"),
    control_samples = c("RAD51_DMSO_26_R1", "RAD51_DMSO_26_R2", "RAD51_DMSO_26_R3"),
    treatment_indices = c(8, 9),
    control_indices = c(22, 23, 24)
  )
)

read_narrowpeak <- function(sample_name) {
  file_path <- file.path(peak_dir, paste0(sample_name, ".macs2_peaks.narrowPeak"))
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  df <- read.delim(file_path, header = FALSE, stringsAsFactors = FALSE)
  GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2] + 1, df[, 3]))
}

read_peak_union <- function(sample_names) {
  peak_list <- lapply(sample_names, read_narrowpeak)
  peak_list <- peak_list[!sapply(peak_list, is.null)]
  if (length(peak_list) == 0) {
    stop("No narrowPeak files were loaded for: ", paste(sample_names, collapse = ", "))
  }
  IRanges::reduce(Reduce(c, peak_list))
}

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

write_significant_bedgraph <- function(res_df, out_dir, prefix) {
  sig_df <- res_df[!is.na(res_df$pvalue) & res_df$pvalue < 0.05, ]
  bedgraph_file <- file.path(out_dir, paste0(prefix, "_significant_log2fc.bedGraph"))
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

run_comparison <- function(definition) {
  cat("\n=== Running comparison:", definition$name, "===\n")

  treatment_union <- read_peak_union(definition$treatment_samples)
  control_union <- read_peak_union(definition$control_samples)

  cat(definition$treatment, "union peaks:", length(treatment_union), "\n")
  cat(definition$control, "union peaks:", length(control_union), "\n")

  final_peaks <- IRanges::reduce(c(treatment_union, control_union))
  cat("Peaks in either group (union):", length(final_peaks), "\n")

  counts_df <- read.delim(counts_file, skip = 1)
  chrom <- counts_df[, 1]
  start <- counts_df[, 2]
  end <- counts_df[, 3]
  merged_gr <- GRanges(seqnames = as.character(chrom), ranges = IRanges(start + 1, end))

  overlaps <- findOverlaps(merged_gr, final_peaks)
  keep_indices <- unique(queryHits(overlaps))
  cat("Peaks from merged table matching consensus:", length(keep_indices), "\n")

  counts_df_filtered <- counts_df[keep_indices, ]
  chrom_filt <- counts_df_filtered[, 1]
  start_filt <- counts_df_filtered[, 2]
  end_filt <- counts_df_filtered[, 3]

  sample_names <- c(definition$treatment_samples, definition$control_samples)
  sample_indices <- c(definition$treatment_indices, definition$control_indices)

  counts_subset <- as.data.frame(apply(counts_df_filtered[, sample_indices], 2, as.integer))
  colnames(counts_subset) <- sample_names
  peak_ids <- paste(chrom_filt, start_filt, end_filt, sep = ":")
  rownames(counts_subset) <- peak_ids
  cat("Final peaks in DESeq2 analysis:", nrow(counts_subset), "\n")

  sample_metadata <- data.frame(
    sample = colnames(counts_subset),
    condition = c(rep(definition$treatment, length(definition$treatment_samples)), rep(definition$control, length(definition$control_samples))),
    row.names = colnames(counts_subset)
  )

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_subset),
    colData = sample_metadata,
    design = ~ condition
  )
  dds$condition <- relevel(dds$condition, ref = definition$control)

  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", definition$treatment, definition$control))

  res_df <- as.data.frame(res)
  res_df$peak_id <- rownames(res_df)
  res_df <- annotate_peaks_with_genes(res_df, chrom_filt, start_filt, end_filt, peak_ids, gtf_file)

  out_dir <- file.path(out_root, definition$name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(res_df, file.path(out_dir, paste0(definition$name, "_results.csv")), row.names = FALSE)

  res_sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
  write.csv(res_sig, file.path(out_dir, paste0(definition$name, "_significant.csv")), row.names = FALSE)
  write_significant_bedgraph(res_df, out_dir, definition$name)

  pdf(file.path(out_dir, "ma_plot.pdf"), width = 8, height = 6)
  plotMA(res, ylim = c(-3, 3))
  dev.off()

  pdf(file.path(out_dir, "volcano_plot.pdf"), width = 8, height = 6)
  volcano_y <- -log10(res_df$pvalue + 1e-300)
  plot(
    res_df$log2FoldChange,
    volcano_y,
    main = paste0("Volcano Plot: ", definition$treatment, " vs ", definition$control),
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

  cat("\n=== DESeq2 Summary ===\n")
  cat("Total peaks analyzed:", nrow(res_df), "\n")
  cat("Significant peaks (padj < 0.05):", nrow(res_sig), "\n")
  cat("Upregulated in", definition$treatment, "(log2FC > 0):", sum(res_df$log2FoldChange > 0 & res_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("Downregulated in", definition$treatment, "(log2FC < 0):", sum(res_df$log2FoldChange < 0 & res_df$padj < 0.05, na.rm = TRUE), "\n")
  cat("Results saved to:", out_dir, "\n")

  invisible(list(out_dir = out_dir, res_df = res_df, res_sig = res_sig))
}

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
results <- lapply(comparison_defs, run_comparison)
