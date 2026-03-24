process CHIPSEEKER_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'
    container 'docker://weishwu/chipseeker:1.42.0'

    input:
    tuple val(meta), path(peaks_annotated)
    path txdb_build
    path genome_annotation

    output:
    tuple val(meta), path("${meta.id}_chipseeker_*"), emit: plots
    tuple val(meta), path("${meta.id}_chipseeker_annotation.tsv"), emit: annotation
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat > chipseeker_annotate.R <<'RSCRIPT'
library(ChIPseeker)
library(ggplot2)

# Create temp library directory for package installation
temp_lib <- file.path(tempdir(), "R_packages")
dir.create(temp_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(temp_lib, .libPaths()))

# Install and load required packages
if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", lib = temp_lib)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}
if (!require("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", lib = temp_lib)
    library(org.Hs.eg.db)
}
if (!require("ggupset", quietly = TRUE)) {
    BiocManager::install("ggupset", lib = temp_lib)
    library(ggupset)
}

# Load annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Plot only canonical chromosomes (exclude scaffolds/random/chrUn)
main_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# Read HOMER annotatePeaks output
peaks_df <- read.delim("${peaks_annotated}", skip=1, stringsAsFactors=FALSE)
colnames(peaks_df) <- c("PeakID", "Chr", "Start", "End", "Strand", "PeakScore", "FocusRatio",
                        "Annotation", "DetailedAnnotation", "DistToTSS", "NearestPromoter",
                        "EntrezID", "UnigeneID", "RefSeqID", "EnsemblID", "GeneName",
                        "GeneAlias", "GeneDescription", "GeneType", "CpG", "GC")

# Convert to GRanges
gr <- GRanges(
    seqnames = peaks_df\$Chr,
    ranges = IRanges(start = as.numeric(peaks_df\$Start), 
                     end = as.numeric(peaks_df\$End)),
    strand = peaks_df\$Strand,
    peakID = peaks_df\$PeakID
)

# Add a score column compatible with covplot(weightCol = "V5")
mcols(gr)\$V5 <- as.numeric(peaks_df\$PeakScore)

# Strictly filter to canonical chromosomes for coverage plotting
gr_cov <- GenomeInfoDb::keepSeqlevels(gr, intersect(GenomeInfoDb::seqlevels(gr), main_chrs), pruning.mode = "coarse")

# Annotate peaks with ChIPseeker
peakAnno <- annotatePeak(
    gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db",
    verbose = FALSE
)

# Save annotation results
anno_result <- as.data.frame(peakAnno)
write.table(anno_result, file="${prefix}_chipseeker_annotation.tsv", 
            sep="\\t", quote=FALSE, row.names=FALSE)

# Create pie chart of genomic annotation
pdf("${prefix}_chipseeker_annotation_pie.pdf", width=8, height=8)
plotAnnoPie(peakAnno)
dev.off()

png("${prefix}_chipseeker_annotation_pie.png", width=800, height=800, res=100)
plotAnnoPie(peakAnno)
dev.off()

# Create bar plot of genomic annotation
pdf("${prefix}_chipseeker_annotation_bar.pdf", width=10, height=6)
plotAnnoBar(peakAnno)
dev.off()

png("${prefix}_chipseeker_annotation_bar.png", width=1000, height=600, res=100)
plotAnnoBar(peakAnno)
dev.off()

# Create upset plot
pdf("${prefix}_chipseeker_annotation_upset.pdf", width=12, height=8)
upsetplot(peakAnno, vennpie = FALSE)
dev.off()

png("${prefix}_chipseeker_annotation_upset.png", width=1200, height=800, res=100)
upsetplot(peakAnno, vennpie = FALSE)
dev.off()

# Distribution of peaks relative to TSS
pdf("${prefix}_chipseeker_dist_to_tss.pdf", width=10, height=6)
plotDistToTSS(peakAnno, title = "Distribution of peaks relative to TSS")
dev.off()

png("${prefix}_chipseeker_dist_to_tss.png", width=1000, height=600, res=100)
plotDistToTSS(peakAnno, title = "Distribution of peaks relative to TSS")
dev.off()

# Coverage of peak regions across chromosomes
pdf("${prefix}_chipseeker_coverage.pdf", width=12, height=6)
covplot(gr_cov, weightCol = "V5", chrs = main_chrs)
dev.off()

png("${prefix}_chipseeker_coverage.png", width=1200, height=600, res=100)
covplot(gr_cov, weightCol = "V5", chrs = main_chrs)
dev.off()

cat("ChIPseeker annotation completed\\n")
RSCRIPT

    Rscript chipseeker_annotate.R

    printf '"ChIPseeker":\n' > versions.yml
    printf '    ChIPseeker: %s\n' '1.42.0' >> versions.yml
    printf '    TxDb.Hsapiens.UCSC.hg38.knownGene: %s\n' '3.20.0' >> versions.yml
    """
}

process CHIPSEEKER_COMPARE {
    label 'process_medium'
    container 'docker://weishwu/chipseeker:1.42.0'

    input:
    path(annotation_files)

    output:
    path("chipseeker_comparison_annotation_bar.pdf"), emit: anno_bar_pdf
    path("chipseeker_comparison_annotation_bar.png"), emit: anno_bar_png
    path("chipseeker_comparison_dist_to_tss.pdf"), emit: dist_pdf
    path("chipseeker_comparison_dist_to_tss.png"), emit: dist_png
    path("chipseeker_comparison_coverage.pdf"), emit: cov_pdf
    path("chipseeker_comparison_coverage.png"), emit: cov_png
    path("chipseeker_comparison_coverage_condition_average.pdf"), emit: cov_avg_pdf
    path("chipseeker_comparison_coverage_condition_average.png"), emit: cov_avg_png
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    cat > chipseeker_compare.R <<'RSCRIPT'
library(ChIPseeker)
library(ggplot2)

# Create temp library directory for package installation
temp_lib <- file.path(tempdir(), "R_packages")
dir.create(temp_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(temp_lib, .libPaths()))

# Install and load required packages
if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", lib = temp_lib)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}
if (!require("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", lib = temp_lib)
    library(org.Hs.eg.db)
}
if (!require("ggupset", quietly = TRUE)) {
    BiocManager::install("ggupset", lib = temp_lib)
    library(ggupset)
}

# Load annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Plot only canonical chromosomes (exclude scaffolds/random/chrUn)
main_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

# List all annotation files
anno_files <- list.files(".", pattern = "_chipseeker_annotation\\\\.tsv\$", full.names = TRUE)

make_placeholder_pdf <- function(path, label) {
    pdf(path, width=12, height=6)
    plot.new()
    text(0.5, 0.5, label)
    dev.off()
}

make_placeholder_png <- function(path, label) {
    png(path, width=1200, height=600, res=100)
    plot.new()
    text(0.5, 0.5, label)
    dev.off()
}

if (length(anno_files) > 0) {
    # Read and collect annotations for comparison
    peakAnnoList <- list()
    anno_all <- list()
    anno_stat_all <- list()
    
    for (f in anno_files) {
        # Read annotation TSV from CHIPSEEKER_ANNOTATE
        anno_df <- read.delim(f, stringsAsFactors = FALSE)
        
        # Convert back to numeric types
        seqnames <- anno_df\$seqnames
        starts <- as.numeric(anno_df\$start)
        ends <- as.numeric(anno_df\$end)
        strands <- anno_df\$strand
        
        # Create GRanges
        gr <- GRanges(
            seqnames = seqnames,
            ranges = IRanges(start = starts, end = ends),
            strand = strands
        )

        # Ensure V5 is present for covplot weighting
        if ("V5" %in% colnames(anno_df)) {
            mcols(gr)\$V5 <- as.numeric(anno_df\$V5)
        } else if ("score" %in% colnames(anno_df)) {
            mcols(gr)\$V5 <- as.numeric(anno_df\$score)
        } else {
            mcols(gr)\$V5 <- rep(1, nrow(anno_df))
        }

        # Strictly filter to canonical chromosomes for coverage plotting
        gr <- GenomeInfoDb::keepSeqlevels(gr, intersect(GenomeInfoDb::seqlevels(gr), main_chrs), pruning.mode = "coarse")

        # Keep flattened annotation table for robust custom comparison plots
        anno_df\$sample <- sub("_chipseeker_annotation\\\\.tsv", "", basename(f))
        anno_all[[basename(f)]] <- anno_df

        # Use ChIPseeker's own annotation pipeline: annotatePeak -> getAnnoStat
        peakAnno_cs <- tryCatch(
            annotatePeak(
                gr,
                tssRegion = c(-3000, 3000),
                TxDb = txdb,
                annoDb = "org.Hs.eg.db",
                verbose = FALSE
            ),
            error = function(e) {
                cat("Note: annotatePeak failed for", anno_df\$sample[1], "-", conditionMessage(e), "\\n")
                NULL
            }
        )

        if (!is.null(peakAnno_cs)) {
            stat_obj <- tryCatch(getAnnoStat(peakAnno_cs), error = function(e) {
                cat("Note: getAnnoStat failed for", anno_df\$sample[1], "-", conditionMessage(e), "\\n")
                NULL
            })
            if (!is.null(stat_obj)) {
                stat_df <- as.data.frame(stat_obj)
                feature_vals <- NULL
                if ("Feature" %in% colnames(stat_df)) {
                    feature_vals <- as.character(stat_df\$Feature)
                } else if (ncol(stat_df) > 0 && !is.numeric(stat_df[[1]])) {
                    feature_vals <- as.character(stat_df[[1]])
                } else if (!is.null(rownames(stat_df)) && !all(grepl("^[0-9]+\$", rownames(stat_df)))) {
                    feature_vals <- rownames(stat_df)
                } else if (!is.null(names(stat_obj)) && length(names(stat_obj)) == nrow(stat_df)) {
                    feature_vals <- names(stat_obj)
                } else {
                    feature_vals <- as.character(seq_len(nrow(stat_df)))
                }

                freq_col <- if ("Frequency" %in% colnames(stat_df)) "Frequency" else {
                    num_cols <- colnames(stat_df)[vapply(stat_df, is.numeric, logical(1))]
                    if (length(num_cols) > 0) num_cols[1] else colnames(stat_df)[1]
                }
                stat_df\$Feature <- feature_vals
                stat_df\$Frequency <- as.numeric(stat_df[[freq_col]])
                stat_df\$sample <- anno_df\$sample[1]
                anno_stat_all[[anno_df\$sample[1]]] <- stat_df[, c("sample", "Feature", "Frequency")]
            }
        }
        
        peakAnnoList[[basename(f)]] <- gr
    }
    
    # Name by sample
    sample_names <- sub("_chipseeker_annotation\\\\.tsv", "", basename(anno_files))
    names(peakAnnoList) <- sample_names
    
    # Comparison plots (only if we have multiple samples)
    if (length(peakAnnoList) > 1) {
        all_df <- do.call(rbind, anno_all)

        # Create comparison annotation bar plot using ChIPseeker getAnnoStat categories
        pdf("chipseeker_comparison_annotation_bar.pdf", width=12, height=6)
        tryCatch({
            if (length(anno_stat_all) == 0) stop("No getAnnoStat output available")
            pdat <- do.call(rbind, anno_stat_all)
            bar_feature_levels <- unique(pdat\$Feature)
            pdat\$Feature <- factor(pdat\$Feature, levels = bar_feature_levels)
            bar_feature_colors <- setNames(ChIPseeker:::getCols(length(bar_feature_levels)), bar_feature_levels)
            p <- ggplot2::ggplot(pdat, ggplot2::aes(x = sample, y = Frequency, fill = Feature)) +
                ggplot2::geom_col(position = "stack") +
                ggplot2::scale_fill_manual(values = bar_feature_colors, breaks = bar_feature_levels) +
                ggplot2::ylab("Frequency (%)") +
                ggplot2::xlab("Sample") +
                ggplot2::ggtitle("Peak Annotation Comparison (ChIPseeker getAnnoStat)") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
            print(p)
        }, error = function(e) {
            cat("Note: Comparison bar plot generation skipped -", conditionMessage(e), "\\n")
        })
        dev.off()
        
        png("chipseeker_comparison_annotation_bar.png", width=1200, height=600, res=100)
        tryCatch({
            if (length(anno_stat_all) == 0) stop("No getAnnoStat output available")
            pdat <- do.call(rbind, anno_stat_all)
            bar_feature_levels <- unique(pdat\$Feature)
            pdat\$Feature <- factor(pdat\$Feature, levels = bar_feature_levels)
            bar_feature_colors <- setNames(ChIPseeker:::getCols(length(bar_feature_levels)), bar_feature_levels)
            p <- ggplot2::ggplot(pdat, ggplot2::aes(x = sample, y = Frequency, fill = Feature)) +
                ggplot2::geom_col(position = "stack") +
                ggplot2::scale_fill_manual(values = bar_feature_colors, breaks = bar_feature_levels) +
                ggplot2::ylab("Frequency (%)") +
                ggplot2::xlab("Sample") +
                ggplot2::ggtitle("Peak Annotation Comparison (ChIPseeker getAnnoStat)") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
            print(p)
        }, error = function(e) {
            cat("Note: Comparison bar plot generation skipped -", conditionMessage(e), "\\n")
        })
        dev.off()
        
        # Distance to TSS comparison using annotation table values
        feature_levels <- c("0 bp-1 kb", "1-3 kb", "3-5 kb", "5-10 kb", "10-100 kb", ">100 kb")
        feature_colors <- c(
            "0 bp-1 kb" = "#9ecae1",
            "1-3 kb" = "#377eb8",
            "3-5 kb" = "#c8ad73",
            "5-10 kb" = "#8cc084",
            "10-100 kb" = "#41bfb5",
            ">100 kb" = "#b784c7"
        )

        pdf("chipseeker_comparison_dist_to_tss.pdf", width=12, height=6)
        tryCatch({
            tss_col <- if ("distanceToTSS" %in% colnames(all_df)) "distanceToTSS" else if ("DistToTSS" %in% colnames(all_df)) "DistToTSS" else NA
            if (is.na(tss_col)) stop("No distanceToTSS column found")
            pdat <- data.frame(sample = all_df\$sample, distanceToTSS = as.numeric(all_df[[tss_col]]), stringsAsFactors = FALSE)
            pdat <- pdat[is.finite(pdat\$distanceToTSS), ]
            pdat\$direction <- ifelse(pdat\$distanceToTSS < 0, "Upstream", "Downstream")
            pdat\$abs_kb <- abs(pdat\$distanceToTSS) / 1000
            pdat\$Feature <- cut(
                pdat\$abs_kb,
                breaks = c(0, 1, 3, 5, 10, 100, Inf),
                labels = c("0 bp-1 kb", "1-3 kb", "3-5 kb", "5-10 kb", "10-100 kb", ">100 kb"),
                right = TRUE,
                include.lowest = TRUE
            )

            sum_df <- as.data.frame(table(pdat\$sample, pdat\$direction, pdat\$Feature), stringsAsFactors = FALSE)
            colnames(sum_df) <- c("sample", "direction", "Feature", "count")
            sum_df <- sum_df[sum_df\$count > 0, ]

            dir_totals <- aggregate(count ~ sample + direction, data = sum_df, FUN = sum)
            sum_df <- merge(sum_df, dir_totals, by = c("sample", "direction"), suffixes = c("", "_dir_total"))
            sum_df\$pct_side <- 50 * sum_df\$count / sum_df\$count_dir_total
            sum_df\$signed_pct <- ifelse(sum_df\$direction == "Upstream", -sum_df\$pct_side, sum_df\$pct_side)
            sum_df\$Feature <- factor(sum_df\$Feature, levels = feature_levels)

            p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = sample, y = signed_pct, fill = Feature)) +
                ggplot2::geom_col(width = 0.85, position = ggplot2::position_stack(reverse = TRUE)) +
                ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
                ggplot2::scale_fill_manual(values = feature_colors, breaks = feature_levels) +
                ggplot2::scale_y_continuous(limits = c(-50, 50), breaks = seq(-50, 50, 10)) +
                ggplot2::ylab("Binding sites (%) (5'->3')") +
                ggplot2::xlab("Sample") +
                ggplot2::ggtitle("Distribution of peaks relative to TSS") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
            print(p)
        }, error = function(e) {
            cat("Note: Comparison distance plot generation skipped\\n")
        })
        dev.off()
        
        png("chipseeker_comparison_dist_to_tss.png", width=1200, height=600, res=100)
        tryCatch({
            tss_col <- if ("distanceToTSS" %in% colnames(all_df)) "distanceToTSS" else if ("DistToTSS" %in% colnames(all_df)) "DistToTSS" else NA
            if (is.na(tss_col)) stop("No distanceToTSS column found")
            pdat <- data.frame(sample = all_df\$sample, distanceToTSS = as.numeric(all_df[[tss_col]]), stringsAsFactors = FALSE)
            pdat <- pdat[is.finite(pdat\$distanceToTSS), ]
            pdat\$direction <- ifelse(pdat\$distanceToTSS < 0, "Upstream", "Downstream")
            pdat\$abs_kb <- abs(pdat\$distanceToTSS) / 1000
            pdat\$Feature <- cut(
                pdat\$abs_kb,
                breaks = c(0, 1, 3, 5, 10, 100, Inf),
                labels = c("0 bp-1 kb", "1-3 kb", "3-5 kb", "5-10 kb", "10-100 kb", ">100 kb"),
                right = TRUE,
                include.lowest = TRUE
            )

            sum_df <- as.data.frame(table(pdat\$sample, pdat\$direction, pdat\$Feature), stringsAsFactors = FALSE)
            colnames(sum_df) <- c("sample", "direction", "Feature", "count")
            sum_df <- sum_df[sum_df\$count > 0, ]

            dir_totals <- aggregate(count ~ sample + direction, data = sum_df, FUN = sum)
            sum_df <- merge(sum_df, dir_totals, by = c("sample", "direction"), suffixes = c("", "_dir_total"))
            sum_df\$pct_side <- 50 * sum_df\$count / sum_df\$count_dir_total
            sum_df\$signed_pct <- ifelse(sum_df\$direction == "Upstream", -sum_df\$pct_side, sum_df\$pct_side)
            sum_df\$Feature <- factor(sum_df\$Feature, levels = feature_levels)

            p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = sample, y = signed_pct, fill = Feature)) +
                ggplot2::geom_col(width = 0.85, position = ggplot2::position_stack(reverse = TRUE)) +
                ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
                ggplot2::scale_fill_manual(values = feature_colors, breaks = feature_levels) +
                ggplot2::scale_y_continuous(limits = c(-50, 50), breaks = seq(-50, 50, 10)) +
                ggplot2::ylab("Binding sites (%) (5'->3')") +
                ggplot2::xlab("Sample") +
                ggplot2::ggtitle("Distribution of peaks relative to TSS") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
            print(p)
        }, error = function(e) {
            cat("Note: Comparison distance plot generation skipped\\n")
        })
        dev.off()

        # Coverage comparison across samples using weighted density per chromosome
        cov_list <- list()
        for (s in names(peakAnnoList)) {
            grs <- peakAnnoList[[s]]
            if (length(grs) == 0) next
            cov_df <- data.frame(
                sample = s,
                chr = as.character(GenomicRanges::seqnames(grs)),
                mid = floor((as.numeric(IRanges::start(grs)) + as.numeric(IRanges::end(grs))) / 2),
                weight = as.numeric(S4Vectors::mcols(grs)\$V5),
                stringsAsFactors = FALSE
            )
            cov_list[[s]] <- cov_df
        }

        cov_all <- do.call(rbind, cov_list)
        if (!is.null(cov_all) && nrow(cov_all) > 0) {
            cov_all <- cov_all[cov_all\$chr %in% main_chrs, ]
            cov_all\$weight[!is.finite(cov_all\$weight) | cov_all\$weight <= 0] <- 1
            cov_all\$condition <- sub("_R[0-9]+$", "", cov_all\$sample)

            p_cov <- ggplot2::ggplot(cov_all, ggplot2::aes(x = mid, color = sample, weight = weight)) +
                ggplot2::geom_density(linewidth = 0.4, adjust = 0.25, na.rm = TRUE) +
                ggplot2::facet_wrap(~ chr, scales = "free", ncol = 4) +
                ggplot2::xlab("Genomic position (bp)") +
                ggplot2::ylab("Weighted density") +
                ggplot2::ggtitle("Peak Density Comparison Across Main Chromosomes") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

            pdf("chipseeker_comparison_coverage.pdf", width=9, height=14)
            print(p_cov)
            dev.off()

            png("chipseeker_comparison_coverage.png", width=1000, height=1800, res=130)
            print(p_cov)
            dev.off()

            # Condition-average coverage using mean per-sample binned signal
            cov_all\$bin <- floor(cov_all\$mid / 1e6) * 1e6
            cov_sample_bin <- aggregate(weight ~ sample + condition + chr + bin, data = cov_all, FUN = sum)
            cov_cond_avg <- aggregate(weight ~ condition + chr + bin, data = cov_sample_bin, FUN = mean)

            p_cov_avg <- ggplot2::ggplot(cov_cond_avg, ggplot2::aes(x = bin, y = weight, color = condition)) +
                ggplot2::geom_line(alpha = 0.95, linewidth = 0.5) +
                ggplot2::facet_wrap(~ chr, scales = "free", ncol = 4) +
                ggplot2::xlab("Genomic position (1 Mb bins)") +
                ggplot2::ylab("Mean weighted peak coverage") +
                ggplot2::ggtitle("Condition-Mean Peak Coverage Across Main Chromosomes") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

            pdf("chipseeker_comparison_coverage_condition_average.pdf", width=9, height=14)
            print(p_cov_avg)
            dev.off()

            png("chipseeker_comparison_coverage_condition_average.png", width=1000, height=1800, res=130)
            print(p_cov_avg)
            dev.off()
        }
    }
}

# Guarantee required comparison files exist for Nextflow output validation
if (!file.exists("chipseeker_comparison_annotation_bar.pdf")) {
    make_placeholder_pdf("chipseeker_comparison_annotation_bar.pdf", "Comparison annotation bar unavailable")
}
if (!file.exists("chipseeker_comparison_annotation_bar.png")) {
    make_placeholder_png("chipseeker_comparison_annotation_bar.png", "Comparison annotation bar unavailable")
}
if (!file.exists("chipseeker_comparison_dist_to_tss.pdf")) {
    make_placeholder_pdf("chipseeker_comparison_dist_to_tss.pdf", "Comparison distance-to-TSS unavailable")
}
if (!file.exists("chipseeker_comparison_dist_to_tss.png")) {
    make_placeholder_png("chipseeker_comparison_dist_to_tss.png", "Comparison distance-to-TSS unavailable")
}
if (!file.exists("chipseeker_comparison_coverage.pdf")) {
    make_placeholder_pdf("chipseeker_comparison_coverage.pdf", "Comparison coverage unavailable")
}
if (!file.exists("chipseeker_comparison_coverage.png")) {
    make_placeholder_png("chipseeker_comparison_coverage.png", "Comparison coverage unavailable")
}
if (!file.exists("chipseeker_comparison_coverage_condition_average.pdf")) {
    make_placeholder_pdf("chipseeker_comparison_coverage_condition_average.pdf", "Condition-average comparison coverage unavailable")
}
if (!file.exists("chipseeker_comparison_coverage_condition_average.png")) {
    make_placeholder_png("chipseeker_comparison_coverage_condition_average.png", "Condition-average comparison coverage unavailable")
}

cat("ChIPseeker comparison completed\\n")
RSCRIPT

    Rscript chipseeker_compare.R

    printf '"ChIPseeker_compare":\n' > versions.yml
    printf '    ChIPseeker: %s\n' '1.42.0' >> versions.yml
    """
}
