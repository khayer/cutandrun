#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(tokens) {
    out <- list()
    key <- NULL
    for (token in tokens) {
        if (startsWith(token, "--")) {
            key <- sub("^--", "", token)
            out[[key]] <- TRUE
        } else if (!is.null(key)) {
            out[[key]] <- token
            key <- NULL
        }
    }
    out
}

opts <- parse_args(args)

required <- c("counts", "annotation", "out-prefix", "group-a", "group-b")
missing <- required[!vapply(required, function(key) !is.null(opts[[key]]), logical(1))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
}

counts_file <- opts[["counts"]]
annotation_file <- opts[["annotation"]]
out_prefix <- opts[["out-prefix"]]
group_a <- opts[["group-a"]]
group_b <- opts[["group-b"]]
fdr_cutoff <- if (!is.null(opts[["fdr"]])) as.numeric(opts[["fdr"]]) else 0.05
log2fc_cutoff <- if (!is.null(opts[["log2fc-cutoff"]])) as.numeric(opts[["log2fc-cutoff"]]) else 1.0
promoter_window <- if (!is.null(opts[["promoter-window"]])) as.integer(opts[["promoter-window"]]) else 3000
top_n <- if (!is.null(opts[["top-n"]])) as.integer(opts[["top-n"]]) else 1000

# Set up writable temp library for package installation in containers
temp_lib <- file.path(tempdir(), "R_packages")
if (!dir.exists(temp_lib)) {
    dir.create(temp_lib, recursive = TRUE)
}
.libPaths(c(temp_lib, .libPaths()))

ensure_package <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = temp_lib, quiet = TRUE)
        }
        if (bioc) {
            BiocManager::install(pkg, ask = FALSE, update = FALSE, lib = temp_lib)
        } else {
            install.packages(pkg, repos = "https://cloud.r-project.org", lib = temp_lib, quiet = TRUE)
        }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

ensure_package("ggplot2")
ensure_package("DESeq2", bioc = TRUE)

prefix <- out_prefix

counts <- read.delim(counts_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)
if (ncol(counts) < 4) {
    stop("Counts table must contain genomic coordinates plus at least one sample column")
}

names(counts)[1:3] <- c("Chr", "Start", "End")
counts$Start <- as.integer(counts$Start)
counts$End <- as.integer(counts$End)

# Strip quotes and # from column names (from deeptools multiBamSummary output)
colnames(counts) <- gsub("^#?'|'$", "", colnames(counts))

annotation <- read.delim(annotation_file, header = FALSE, skip = 1, sep = "\t", check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
expected_names <- c(
    "PeakID", "Chr", "Start", "End", "Strand", "PeakScore", "FocusRatio", "Annotation",
    "DetailedAnnotation", "DistanceToTSS", "NearestPromoterID", "EntrezID", "NearestUnigene",
    "NearestRefseq", "NearestEnsembl", "GeneName", "GeneAlias", "GeneDescription", "GeneType",
    "CpG_percent", "GC_percent"
)
if (ncol(annotation) < 21) {
    stop(sprintf("Unexpected annotatePeaks table format: found %d columns, expected at least 21", ncol(annotation)))
}
cat(sprintf("Annotation table has %d columns, assigning names to first %d\n", ncol(annotation), length(expected_names)))
colnames(annotation) <- expected_names[1:ncol(annotation)]
annotation$Start <- as.integer(annotation$Start)
annotation$End <- as.integer(annotation$End)
annotation$DistanceToTSS <- suppressWarnings(as.numeric(annotation$DistanceToTSS))
annotation$GC_percent <- suppressWarnings(as.numeric(annotation$GC_percent))

# Counts are from BED format (0-based), annotation is from HOMER (1-based)
# Convert counts to 1-based to match annotation
counts$Start <- as.integer(counts$Start) + 1

counts_cols <- setdiff(names(counts), c("Chr", "Start", "End"))
group_a_cols <- counts_cols[startsWith(counts_cols, group_a)]
group_b_cols <- counts_cols[startsWith(counts_cols, group_b)]
selected_cols <- c(group_a_cols, group_b_cols)

cat("\nGroup matching debug:\n")
cat(sprintf("  All sample columns: %s\n", paste(head(counts_cols, 10), collapse=", ")))
cat(sprintf("  Group A ('%s') columns found: %d\n", group_a, length(group_a_cols)))
if (length(group_a_cols) > 0) {
    cat(sprintf("    %s\n", paste(group_a_cols, collapse=", ")))
}
cat(sprintf("  Group B ('%s') columns found: %d\n", group_b, length(group_b_cols)))
if (length(group_b_cols) > 0) {
    cat(sprintf("    %s\n", paste(group_b_cols, collapse=", ")))
}

if (length(group_a_cols) < 2) {
    stop(sprintf("Need at least two samples for group_a '%s'", group_a))
}
if (length(group_b_cols) < 2) {
    stop(sprintf("Need at least two samples for group_b '%s'", group_b))
}

counts <- counts[, c("Chr", "Start", "End", selected_cols), drop = FALSE]

peak_key <- function(df) paste(df$Chr, df$Start, df$End, sep = ":")
count_key <- peak_key(counts)
annot_key <- peak_key(annotation)

# Merge annotation and counts tables on genomic coordinates
# Add keys for merging
counts$peak_id <- count_key
annotation$peak_id <- annot_key

# Debug: Show first few keys and data types
cat("Annotation data after processing:\n")
cat(sprintf("  Rows: %d, Cols: %d\n", nrow(annotation), ncol(annotation)))
cat(sprintf("  First peak: %s-%s, Start class: %s, End class: %s\n", 
            annotation$Chr[1], annotation$PeakID[1], class(annotation$Start[1]), class(annotation$End[1])))
cat(sprintf("  First peak key: %s\n", annot_key[1]))

cat("\nCounts data after processing:\n")
cat(sprintf("  Rows: %d, Cols (first few): %s\n", nrow(counts), paste(colnames(counts)[1:5], collapse=", ")))
cat(sprintf("  First row Chr: %s, Start: %s, End: %s (classes: %s, %s, %s)\n",
            counts$Chr[1], counts$Start[1], counts$End[1],
            class(counts$Chr[1]), class(counts$Start[1]), class(counts$End[1])))
cat(sprintf("  First row key: %s\n", count_key[1]))

cat("\nMerge keys - checking first 3 from each:\n")
cat("Annotation keys:\n")
print(head(annot_key, 3))
cat("Counts keys:\n")
print(head(count_key, 3))

# Merge on peak_id
annotated <- merge(annotation, counts[, c("peak_id", selected_cols), drop = FALSE], by = "peak_id", all = FALSE)

if (nrow(annotated) == 0) {
    # More detailed error message
    n_ann_unique <- length(unique(annot_key))
    n_count_unique <- length(unique(count_key))
    n_ann_in_counts <- sum(annotation$peak_id %in% counts$peak_id)
    cat(sprintf("Annotations: %d total, %d unique peaks\n", nrow(annotation), n_ann_unique))
    cat(sprintf("Counts: %d total, %d unique peaks\n", nrow(counts), n_count_unique))
    cat(sprintf("Peaks from annotation found in counts: %d\n", n_ann_in_counts))
    stop("No matching peaks found between counts and annotation tables")
}

if (nrow(annotated) < nrow(annotation)) {
    warning(sprintf("Only %d of %d annotation peaks found in counts table", nrow(annotated), nrow(annotation)))
}

count_matrix <- as.matrix(annotated[, selected_cols, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
count_matrix <- round(count_matrix)
rownames(count_matrix) <- annotated$peak_id

condition <- ifelse(startsWith(selected_cols, group_a), group_a, group_b)
sample_table <- data.frame(row.names = selected_cols, condition = factor(condition, levels = c(group_b, group_a)))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_table, design = ~ condition)

# Before filtering, get the peak_ids that have zero counts
peaks_with_counts <- names(which(rowSums(DESeq2::counts(dds)) > 0))
peaks_to_remove <- names(which(rowSums(DESeq2::counts(dds)) == 0))

cat(sprintf("Total peaks in DDS: %d\n", nrow(dds)))
cat(sprintf("Peaks with zero counts: %d\n", length(peaks_to_remove)))
cat(sprintf("Peaks to keep: %d\n", length(peaks_with_counts)))

# Subset annotated dataframe to keep only peaks with counts
# Match by peak_id to ensure correct correspondence
annotated_kept <- annotated[annotated$peak_id %in% peaks_with_counts, , drop = FALSE]

cat(sprintf("Filtering out %d peaks with zero counts, keeping %d peaks\n", 
            nrow(annotated) - nrow(annotated_kept), nrow(annotated_kept)))

# Filter DDS to keep only peaks with counts
dds <- dds[peaks_with_counts, ]
dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds, contrast = c("condition", group_a, group_b), alpha = fdr_cutoff)
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

cat(sprintf("DESeq2 results: %d peaks\n", nrow(res_df)))
cat(sprintf("Annotated peaks kept: %d peaks\n", nrow(annotated_kept)))

result_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

# Initialize result columns in the annotated_kept dataframe
for (col in result_cols) {
    annotated_kept[[col]] <- NA_real_
}

# Match DESeq2 results to annotated peaks
# Use annotated_kept for matching since those are the ones DESeq2 saw
match_index <- match(annotated_kept$peak_id, res_df$peak_id)
if (all(is.na(match_index))) {
    stop(sprintf("Could not match any DESeq2 results to annotated peaks. DESeq2 had %d peaks with rownames, annotated_kept had %d peaks with peak_id", nrow(res_df), nrow(annotated_kept)))
}
present <- !is.na(match_index)

# Add results to the annotated_kept dataframe
for (col in result_cols) {
    annotated_kept[[col]][present] <- res_df[[col]][match_index[present]]
}

annotated_kept$promoter_window <- !is.na(annotated_kept$DistanceToTSS) & abs(annotated_kept$DistanceToTSS) <= promoter_window
annotated_kept$significant <- !is.na(annotated_kept$padj) & annotated_kept$padj <= fdr_cutoff
annotated_kept$status <- "non_promoter"
annotated_kept$status[annotated_kept$promoter_window & !annotated_kept$significant] <- "unchanged"
annotated_kept$status[annotated_kept$promoter_window & annotated_kept$significant & annotated_kept$log2FoldChange <= -log2fc_cutoff] <- "loss"
annotated_kept$status[annotated_kept$promoter_window & annotated_kept$significant & annotated_kept$log2FoldChange >= log2fc_cutoff] <- "gain"
annotated_kept$status <- factor(annotated_kept$status, levels = c("loss", "unchanged", "gain", "non_promoter"))

# Write results for all peaks that passed DESeq2 filtering (had non-zero counts)
write.table(annotated_kept, file = paste0(prefix, ".deseq2_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

promoter_df <- annotated_kept[annotated_kept$promoter_window, , drop = FALSE]
write.table(promoter_df, file = paste0(prefix, ".promoter_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

promoter_loss <- promoter_df[promoter_df$status == "loss", , drop = FALSE]
promoter_unchanged <- promoter_df[promoter_df$status == "unchanged", , drop = FALSE]
write.table(promoter_loss, file = paste0(prefix, ".promoter_loss.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(promoter_unchanged, file = paste0(prefix, ".promoter_not_affected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

gc_df <- promoter_df[promoter_df$status %in% c("loss", "unchanged") & !is.na(promoter_df$GC_percent), , drop = FALSE]

gc_groups <- split(gc_df, droplevels(gc_df$status))
gc_groups <- gc_groups[lengths(gc_groups) > 0]
gc_summary <- if (length(gc_groups) == 0) {
    gc_summary <- data.frame(status = character(), n_peaks = integer(), mean_gc_percent = numeric(), median_gc_percent = numeric(), sd_gc_percent = numeric())
} else {
    do.call(rbind, lapply(names(gc_groups), function(status_name) {
        df <- gc_groups[[status_name]]
        data.frame(
            status = status_name,
            n_peaks = nrow(df),
            mean_gc_percent = mean(df$GC_percent),
            median_gc_percent = median(df$GC_percent),
            sd_gc_percent = if (nrow(df) > 1) stats::sd(df$GC_percent) else NA_real_,
            stringsAsFactors = FALSE
        )
    }))
}
if (is.null(gc_summary) || nrow(gc_summary) == 0) {
    gc_summary <- data.frame(status = character(), n_peaks = integer(), mean_gc_percent = numeric(), median_gc_percent = numeric(), sd_gc_percent = numeric())
}
write.table(gc_summary, file = paste0(prefix, ".promoter_gc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (nrow(gc_df) >= 2 && length(unique(gc_df$status)) == 2) {
    gc_test <- stats::wilcox.test(GC_percent ~ status, data = gc_df, exact = FALSE)
    gc_test_df <- data.frame(
        comparison = paste(group_a, "vs", group_b),
        promoter_window_bp = promoter_window,
        fdr_cutoff = fdr_cutoff,
        log2fc_cutoff = log2fc_cutoff,
        n_loss = sum(gc_df$status == "loss"),
        n_unchanged = sum(gc_df$status == "unchanged"),
        statistic = unname(gc_test$statistic),
        p_value = gc_test$p.value,
        median_gc_loss = median(gc_df$GC_percent[gc_df$status == "loss"]),
        median_gc_unchanged = median(gc_df$GC_percent[gc_df$status == "unchanged"]),
        median_difference = median(gc_df$GC_percent[gc_df$status == "loss"]) - median(gc_df$GC_percent[gc_df$status == "unchanged"]),
        stringsAsFactors = FALSE
    )
} else {
    gc_test_df <- data.frame(
        comparison = paste(group_a, "vs", group_b),
        promoter_window_bp = promoter_window,
        fdr_cutoff = fdr_cutoff,
        log2fc_cutoff = log2fc_cutoff,
        n_loss = sum(gc_df$status == "loss"),
        n_unchanged = sum(gc_df$status == "unchanged"),
        statistic = NA_real_,
        p_value = NA_real_,
        median_gc_loss = if (sum(gc_df$status == "loss") > 0) median(gc_df$GC_percent[gc_df$status == "loss"]) else NA_real_,
        median_gc_unchanged = if (sum(gc_df$status == "unchanged") > 0) median(gc_df$GC_percent[gc_df$status == "unchanged"]) else NA_real_,
        median_difference = if (sum(gc_df$status == "loss") > 0 && sum(gc_df$status == "unchanged") > 0) {
            median(gc_df$GC_percent[gc_df$status == "loss"]) - median(gc_df$GC_percent[gc_df$status == "unchanged"])
        } else {
            NA_real_
        },
        stringsAsFactors = FALSE
    )
}
write.table(gc_test_df, file = paste0(prefix, ".promoter_gc_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# Additional ranking-based outputs using raw p-values.
# This is useful when the adjusted p-value threshold is too strict for small replicate numbers.
ranked_promoter_df <- promoter_df[!is.na(promoter_df$pvalue) & !is.na(promoter_df$log2FoldChange), , drop = FALSE]
ranked_loss <- ranked_promoter_df[ranked_promoter_df$log2FoldChange < 0, , drop = FALSE]
ranked_gain <- ranked_promoter_df[ranked_promoter_df$log2FoldChange > 0, , drop = FALSE]

ranked_loss <- ranked_loss[order(ranked_loss$pvalue, -abs(ranked_loss$log2FoldChange), -ranked_loss$baseMean), , drop = FALSE]
ranked_gain <- ranked_gain[order(ranked_gain$pvalue, -abs(ranked_gain$log2FoldChange), -ranked_gain$baseMean), , drop = FALSE]

if (nrow(ranked_loss) > top_n) {
    ranked_loss <- ranked_loss[seq_len(top_n), , drop = FALSE]
}
if (nrow(ranked_gain) > top_n) {
    ranked_gain <- ranked_gain[seq_len(top_n), , drop = FALSE]
}

write.table(ranked_loss, file = paste0(prefix, ".top", top_n, "_loss.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(ranked_gain, file = paste0(prefix, ".top", top_n, "_gain.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

ranked_gc_df <- rbind(
    transform(ranked_loss, status = "loss"),
    transform(ranked_gain, status = "gain")
)

if (nrow(ranked_gc_df) == 0) {
    ranked_gc_summary <- data.frame(status = character(), n_peaks = integer(), mean_gc_percent = numeric(), median_gc_percent = numeric(), sd_gc_percent = numeric())
    ranked_gc_test_df <- data.frame(
        comparison = paste0("top", top_n, " raw p-value ranked loss vs gain"),
        top_n = top_n,
        n_loss = nrow(ranked_loss),
        n_gain = nrow(ranked_gain),
        statistic = NA_real_,
        p_value = NA_real_,
        median_gc_loss = if (nrow(ranked_loss) > 0) median(ranked_loss$GC_percent, na.rm = TRUE) else NA_real_,
        median_gc_gain = if (nrow(ranked_gain) > 0) median(ranked_gain$GC_percent, na.rm = TRUE) else NA_real_,
        median_difference = if (nrow(ranked_loss) > 0 && nrow(ranked_gain) > 0) {
            median(ranked_loss$GC_percent, na.rm = TRUE) - median(ranked_gain$GC_percent, na.rm = TRUE)
        } else {
            NA_real_
        },
        stringsAsFactors = FALSE
    )
} else {
    ranked_gc_summary <- do.call(rbind, lapply(split(ranked_gc_df, ranked_gc_df$status), function(df) {
        data.frame(
            status = as.character(unique(df$status)),
            n_peaks = nrow(df),
            mean_gc_percent = mean(df$GC_percent),
            median_gc_percent = median(df$GC_percent),
            sd_gc_percent = if (nrow(df) > 1) stats::sd(df$GC_percent) else NA_real_,
            stringsAsFactors = FALSE
        )
    }))

    if (nrow(ranked_loss) >= 2 && nrow(ranked_gain) >= 2) {
        ranked_gc_test <- stats::wilcox.test(ranked_loss$GC_percent, ranked_gain$GC_percent, exact = FALSE)
        ranked_gc_test_df <- data.frame(
            comparison = paste0("top", top_n, " raw p-value ranked loss vs gain"),
            top_n = top_n,
            n_loss = nrow(ranked_loss),
            n_gain = nrow(ranked_gain),
            statistic = unname(ranked_gc_test$statistic),
            p_value = ranked_gc_test$p.value,
            median_gc_loss = median(ranked_loss$GC_percent),
            median_gc_gain = median(ranked_gain$GC_percent),
            median_difference = median(ranked_loss$GC_percent) - median(ranked_gain$GC_percent),
            stringsAsFactors = FALSE
        )
    } else {
        ranked_gc_test_df <- data.frame(
            comparison = paste0("top", top_n, " raw p-value ranked loss vs gain"),
            top_n = top_n,
            n_loss = nrow(ranked_loss),
            n_gain = nrow(ranked_gain),
            statistic = NA_real_,
            p_value = NA_real_,
            median_gc_loss = if (nrow(ranked_loss) > 0) median(ranked_loss$GC_percent) else NA_real_,
            median_gc_gain = if (nrow(ranked_gain) > 0) median(ranked_gain$GC_percent) else NA_real_,
            median_difference = if (nrow(ranked_loss) > 0 && nrow(ranked_gain) > 0) {
                median(ranked_loss$GC_percent) - median(ranked_gain$GC_percent)
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }
}

if (is.null(ranked_gc_summary) || nrow(ranked_gc_summary) == 0) {
    ranked_gc_summary <- data.frame(status = character(), n_peaks = integer(), mean_gc_percent = numeric(), median_gc_percent = numeric(), sd_gc_percent = numeric())
}

write.table(ranked_gc_summary, file = paste0(prefix, ".top", top_n, "_promoter_gc_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(ranked_gc_test_df, file = paste0(prefix, ".top", top_n, "_promoter_gc_test.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_gc <- function(df, title_suffix) {
    if (nrow(df) == 0) {
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No promoter GC data available"))
    }
    ggplot2::ggplot(df, ggplot2::aes(x = status, y = GC_percent, fill = status)) +
        ggplot2::geom_boxplot(width = 0.55, outlier.alpha = 0.35) +
        ggplot2::geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
        ggplot2::scale_fill_manual(values = c(loss = "#D95F02", unchanged = "#1B9E77", gain = "#7570B3", non_promoter = "#999999"), drop = FALSE) +
        ggplot2::labs(
            title = paste0("Promoter GC content: ", group_a, " vs ", group_b),
            subtitle = title_suffix,
            x = "Promoter status",
            y = "GC content (%)"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
}

plot_gc_violin <- function(df, title_suffix) {
    if (nrow(df) == 0) {
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No promoter GC data available"))
    }
    ggplot2::ggplot(df, ggplot2::aes(x = status, y = GC_percent, fill = status)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.65) +
        ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.85) +
        ggplot2::scale_fill_manual(values = c(loss = "#D95F02", unchanged = "#1B9E77", gain = "#7570B3", non_promoter = "#999999"), drop = FALSE) +
        ggplot2::labs(
            title = paste0("Promoter GC content (violin): ", group_a, " vs ", group_b),
            subtitle = title_suffix,
            x = "Promoter status",
            y = "GC content (%)"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
}

plot_volcano <- function(df) {
    df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
    ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)) +
        ggplot2::geom_point(alpha = 0.45, size = 1, aes(color = status)) +
        ggplot2::scale_color_manual(values = c(loss = "#D95F02", unchanged = "#1B9E77", gain = "#7570B3", non_promoter = "#999999"), drop = FALSE) +
        ggplot2::labs(
            title = paste0("Differential binding volcano: ", group_a, " vs ", group_b),
            x = "log2 fold change",
            y = "-log10 adjusted p-value"
        ) +
        ggplot2::theme_bw()
}

plot_ma <- function(df) {
    ggplot2::ggplot(df, ggplot2::aes(x = log10(baseMean + 1), y = log2FoldChange)) +
        ggplot2::geom_point(alpha = 0.45, size = 1, aes(color = status)) +
        ggplot2::scale_color_manual(values = c(loss = "#D95F02", unchanged = "#1B9E77", gain = "#7570B3", non_promoter = "#999999"), drop = FALSE) +
        ggplot2::labs(
            title = paste0("Differential binding MA plot: ", group_a, " vs ", group_b),
            x = "log10(mean normalized count + 1)",
            y = "log2 fold change"
        ) +
        ggplot2::theme_bw()
}

ggplot2::ggsave(paste0(prefix, ".promoter_gc_boxplot.png"), plot_gc(gc_df, paste0("Promoter window: ", promoter_window, " bp; FDR <= ", fdr_cutoff)), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".promoter_gc_boxplot.pdf"), plot_gc(gc_df, paste0("Promoter window: ", promoter_window, " bp; FDR <= ", fdr_cutoff)), width = 7.5, height = 5.5)
ggplot2::ggsave(paste0(prefix, ".promoter_gc_violin.png"), plot_gc_violin(gc_df, paste0("Promoter window: ", promoter_window, " bp; FDR <= ", fdr_cutoff)), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".promoter_gc_violin.pdf"), plot_gc_violin(gc_df, paste0("Promoter window: ", promoter_window, " bp; FDR <= ", fdr_cutoff)), width = 7.5, height = 5.5)
ggplot2::ggsave(paste0(prefix, ".volcano.png"), plot_volcano(annotated_kept), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".volcano.pdf"), plot_volcano(annotated_kept), width = 7.5, height = 5.5)
ggplot2::ggsave(paste0(prefix, ".ma_plot.png"), plot_ma(annotated_kept), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".ma_plot.pdf"), plot_ma(annotated_kept), width = 7.5, height = 5.5)

ggplot2::ggsave(paste0(prefix, ".top", top_n, "_promoter_gc_boxplot.png"), plot_gc(ranked_gc_df, paste0("Top ", top_n, " promoter peaks ranked by raw p-value")), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".top", top_n, "_promoter_gc_boxplot.pdf"), plot_gc(ranked_gc_df, paste0("Top ", top_n, " promoter peaks ranked by raw p-value")), width = 7.5, height = 5.5)
ggplot2::ggsave(paste0(prefix, ".top", top_n, "_promoter_gc_violin.png"), plot_gc_violin(ranked_gc_df, paste0("Top ", top_n, " promoter peaks ranked by raw p-value")), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".top", top_n, "_promoter_gc_violin.pdf"), plot_gc_violin(ranked_gc_df, paste0("Top ", top_n, " promoter peaks ranked by raw p-value")), width = 7.5, height = 5.5)

versions <- c(
    '"PROMOTER_GC_DIFFBIND":',
    sprintf('  R: "%s"', getRversion()),
    sprintf('  DESeq2: "%s"', as.character(packageVersion("DESeq2"))),
    sprintf('  ggplot2: "%s"', as.character(packageVersion("ggplot2")))
)
writeLines(versions, con = "versions.yml")