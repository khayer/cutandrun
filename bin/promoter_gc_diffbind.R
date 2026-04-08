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

ensure_package <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        if (bioc) {
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
            install.packages(pkg, repos = "https://cloud.r-project.org")
        }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

ensure_package("ggplot2")
ensure_package("DESeq2", bioc = TRUE)

prefix <- out_prefix
annotation_copy <- paste0(prefix, ".merged_peaks.annotatePeaks.txt")
file.copy(annotation_file, annotation_copy, overwrite = TRUE)

counts <- read.delim(counts_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
if (ncol(counts) < 4) {
    stop("Counts table must contain genomic coordinates plus at least one sample column")
}

names(counts)[1:3] <- c("Chr", "Start", "End")
counts$Start <- as.integer(counts$Start)
counts$End <- as.integer(counts$End)

annotation <- read.delim(annotation_copy, header = FALSE, skip = 1, sep = "\t", check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
expected_names <- c(
    "PeakID", "Chr", "Start", "End", "Strand", "PeakScore", "FocusRatio", "Annotation",
    "DetailedAnnotation", "DistanceToTSS", "NearestPromoterID", "EntrezID", "NearestUnigene",
    "NearestRefseq", "NearestEnsembl", "GeneName", "GeneAlias", "GeneDescription", "GeneType",
    "CpG_percent", "GC_percent"
)
if (ncol(annotation) < 21) {
    stop(sprintf("Unexpected annotatePeaks table format: found %d columns", ncol(annotation)))
}
colnames(annotation)[seq_along(expected_names)] <- expected_names
annotation$Start <- as.integer(annotation$Start)
annotation$End <- as.integer(annotation$End)
annotation$DistanceToTSS <- suppressWarnings(as.numeric(annotation$DistanceToTSS))
annotation$GC_percent <- suppressWarnings(as.numeric(annotation$GC_percent))

counts_cols <- setdiff(names(counts), c("Chr", "Start", "End"))
group_a_cols <- counts_cols[startsWith(counts_cols, group_a)]
group_b_cols <- counts_cols[startsWith(counts_cols, group_b)]
selected_cols <- c(group_a_cols, group_b_cols)

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

if (!identical(count_key, annot_key)) {
    ann_index <- match(count_key, annot_key)
    if (any(is.na(ann_index))) {
        stop("Merged peak coordinates in counts and annotation tables do not align")
    }
    annotation <- annotation[ann_index, , drop = FALSE]
}

annotated <- cbind(annotation, counts[, selected_cols, drop = FALSE])
annotated$peak_id <- peak_key(annotated)

count_matrix <- as.matrix(annotated[, selected_cols, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
count_matrix <- round(count_matrix)
rownames(count_matrix) <- annotated$peak_id

condition <- ifelse(startsWith(selected_cols, group_a), group_a, group_b)
sample_table <- data.frame(row.names = selected_cols, condition = factor(condition, levels = c(group_b, group_a)))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_table, design = ~ condition)
dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]
dds <- DESeq2::DESeq(dds)

res <- DESeq2::results(dds, contrast = c("condition", group_a, group_b), alpha = fdr_cutoff)
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)

result_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
for (col in result_cols) {
    annotated[[col]] <- NA_real_
}

match_index <- match(annotated$peak_id, res_df$peak_id)
present <- !is.na(match_index)
for (col in result_cols) {
    annotated[[col]][present] <- res_df[[col]][match_index[present]]
}

annotated$promoter_window <- !is.na(annotated$DistanceToTSS) & abs(annotated$DistanceToTSS) <= promoter_window
annotated$significant <- !is.na(annotated$padj) & annotated$padj <= fdr_cutoff
annotated$status <- "non_promoter"
annotated$status[annotated$promoter_window & !annotated$significant] <- "unchanged"
annotated$status[annotated$promoter_window & annotated$significant & annotated$log2FoldChange <= -log2fc_cutoff] <- "loss"
annotated$status[annotated$promoter_window & annotated$significant & annotated$log2FoldChange >= log2fc_cutoff] <- "gain"
annotated$status <- factor(annotated$status, levels = c("loss", "unchanged", "gain", "non_promoter"))

write.table(annotated, file = paste0(prefix, ".deseq2_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

promoter_df <- annotated[annotated$promoter_window, , drop = FALSE]
write.table(promoter_df, file = paste0(prefix, ".promoter_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

promoter_loss <- promoter_df[promoter_df$status == "loss", , drop = FALSE]
promoter_unchanged <- promoter_df[promoter_df$status == "unchanged", , drop = FALSE]
write.table(promoter_loss, file = paste0(prefix, ".promoter_loss.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(promoter_unchanged, file = paste0(prefix, ".promoter_not_affected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

gc_df <- promoter_df[promoter_df$status %in% c("loss", "unchanged") & !is.na(promoter_df$GC_percent), , drop = FALSE]

gc_summary <- do.call(rbind, lapply(split(gc_df, gc_df$status), function(df) {
    data.frame(
        status = as.character(unique(df$status)),
        n_peaks = nrow(df),
        mean_gc_percent = mean(df$GC_percent),
        median_gc_percent = median(df$GC_percent),
        sd_gc_percent = if (nrow(df) > 1) stats::sd(df$GC_percent) else NA_real_,
        stringsAsFactors = FALSE
    )
}))
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

plot_gc <- function(df, title_suffix) {
    if (nrow(df) == 0) {
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No promoter GC data available"))
    }
    ggplot2::ggplot(df, ggplot2::aes(x = status, y = GC_percent, fill = status)) +
        ggplot2::geom_boxplot(width = 0.55, outlier.alpha = 0.35) +
        ggplot2::geom_jitter(width = 0.12, alpha = 0.45, size = 1) +
        ggplot2::scale_fill_manual(values = c(loss = "#D95F02", unchanged = "#1B9E77"), drop = FALSE) +
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
        ggplot2::scale_fill_manual(values = c(loss = "#D95F02", unchanged = "#1B9E77"), drop = FALSE) +
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
ggplot2::ggsave(paste0(prefix, ".volcano.png"), plot_volcano(annotated), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".volcano.pdf"), plot_volcano(annotated), width = 7.5, height = 5.5)
ggplot2::ggsave(paste0(prefix, ".ma_plot.png"), plot_ma(annotated), width = 7.5, height = 5.5, dpi = 200)
ggplot2::ggsave(paste0(prefix, ".ma_plot.pdf"), plot_ma(annotated), width = 7.5, height = 5.5)

versions <- c(
    '"PROMOTER_GC_DIFFBIND":',
    sprintf('  R: "%s"', getRversion()),
    sprintf('  DESeq2: "%s"', as.character(packageVersion("DESeq2"))),
    sprintf('  ggplot2: "%s"', as.character(packageVersion("ggplot2")))
)
writeLines(versions, con = "versions.yml")