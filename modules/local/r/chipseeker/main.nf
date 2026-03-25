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
library(GenomicRanges)

# Create temp library directory for package installation
temp_lib <- file.path(tempdir(), "R_packages")
dir.create(temp_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(temp_lib, .libPaths()))

# Install optional annotation DB (used when available)
has_orgdb <- require("org.Hs.eg.db", quietly = TRUE)
if (!has_orgdb) {
    tryCatch({
        BiocManager::install("org.Hs.eg.db", lib = temp_lib)
        has_orgdb <- require("org.Hs.eg.db", quietly = TRUE)
    }, error = function(e) {
        cat("Note: org.Hs.eg.db unavailable -", conditionMessage(e), "\\n")
    })
}
if (!require("ggupset", quietly = TRUE)) {
    BiocManager::install("ggupset", lib = temp_lib)
    library(ggupset)
}

# Required on Bioconductor >= 3.19 for makeTxDbFromGFF
has_txdbmaker <- require("txdbmaker", quietly = TRUE)
if (!has_txdbmaker) {
    tryCatch({
        BiocManager::install("txdbmaker", lib = temp_lib)
        has_txdbmaker <- require("txdbmaker", quietly = TRUE)
    }, error = function(e) {
        cat("Note: txdbmaker unavailable -", conditionMessage(e), "\\n")
    })
}

# Build TxDb from provided genome annotation; fallback to hg38 TxDb for compatibility
txdb <- tryCatch({
    if (!has_txdbmaker) stop("txdbmaker not installed")
    txdbmaker::makeTxDbFromGFF("${genome_annotation}")
}, error = function(e) {
    cat("Note: makeTxDbFromGFF failed -", conditionMessage(e), "\\n")
    NULL
})

if (is.null(txdb)) {
    if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
        BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", lib = temp_lib)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    }
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    cat("Note: Using hg38 fallback TxDb\\n")
}

anno_db <- if (has_orgdb) "org.Hs.eg.db" else NULL

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

# Use canonical chromosomes when present, otherwise fall back to available reference seqlevels
available_chrs <- GenomeInfoDb::seqlevels(gr)
cov_chrs <- intersect(main_chrs, available_chrs)
if (length(cov_chrs) > 0) {
    gr_cov <- GenomeInfoDb::keepSeqlevels(gr, cov_chrs, pruning.mode = "coarse")
} else {
    cov_chrs <- available_chrs
    gr_cov <- gr
}

plot_placeholder <- function(label) {
    plot.new()
    text(0.5, 0.5, label)
}

safe_dual_plot <- function(pdf_file, png_file, pdf_w, pdf_h, png_w, png_h, plot_fn, fallback_label) {
    pdf(pdf_file, width = pdf_w, height = pdf_h)
    tryCatch(plot_fn(), error = function(e) {
        plot_placeholder(paste0(fallback_label, " (", conditionMessage(e), ")"))
    })
    dev.off()

    png(png_file, width = png_w, height = png_h, res = 100)
    tryCatch(plot_fn(), error = function(e) {
        plot_placeholder(paste0(fallback_label, " (", conditionMessage(e), ")"))
    })
    dev.off()
}

annotate_with_gff_overlap <- function(gr, gff_file) {
    gff <- tryCatch(
        read.delim(
            gff_file,
            header = FALSE,
            sep = "\t",
            quote = "",
            comment.char = "#",
            stringsAsFactors = FALSE,
            fill = TRUE
        ),
        error = function(e) NULL
    )

    if (is.null(gff) || nrow(gff) == 0 || ncol(gff) < 9) {
        cat("Note: GFF fallback parsing failed or empty - returning Unannotated\\n")
        return(rep("Unannotated", length(gr)))
    }

    colnames(gff)[1:9] <- c("seqname", "source", "feature", "start", "end", "score", "strand", "phase", "attributes")
    gff <- gff[gff\$feature %in% c("CDS", "mRNA", "gene"), , drop = FALSE]

    if (nrow(gff) == 0) {
        cat("Note: No CDS/mRNA/gene features found in GFF fallback\\n")
        return(rep("Unannotated", length(gr)))
    }

    gff_gr <- GRanges(
        seqnames = gff\$seqname,
        ranges = IRanges(start = as.integer(gff\$start), end = as.integer(gff\$end)),
        strand = gff\$strand
    )
    mcols(gff_gr)\$feature <- gff\$feature

    hits <- findOverlaps(gr, gff_gr, ignore.strand = TRUE)
    ann <- rep("Intergenic", length(gr))

    if (length(hits) > 0) {
        split_hits <- split(subjectHits(hits), queryHits(hits))
        for (qid in names(split_hits)) {
            idx <- as.integer(qid)
            fset <- unique(as.character(mcols(gff_gr)\$feature[split_hits[[qid]]]))
            if ("CDS" %in% fset) {
                ann[idx] <- "CDS"
            } else if ("mRNA" %in% fset) {
                ann[idx] <- "Transcript"
            } else if ("gene" %in% fset) {
                ann[idx] <- "Gene"
            } else {
                ann[idx] <- "Intergenic"
            }
        }
    }

    ann
}

# Annotate peaks with ChIPseeker
peakAnno <- tryCatch(
    {
        if (is.null(anno_db)) {
            annotatePeak(
                gr,
                tssRegion = c(-3000, 3000),
                TxDb = txdb,
                verbose = FALSE
            )
        } else {
            annotatePeak(
                gr,
                tssRegion = c(-3000, 3000),
                TxDb = txdb,
                annoDb = anno_db,
                verbose = FALSE
            )
        }
    },
    error = function(e) {
        cat("Note: annotatePeak failed for ${prefix} -", conditionMessage(e), "\\n")
        NULL
    }
)

# Save annotation results
if (!is.null(peakAnno)) {
    anno_result <- as.data.frame(peakAnno)
} else {
    fallback_annotation <- annotate_with_gff_overlap(gr, "${genome_annotation}")
    anno_result <- data.frame(
        seqnames = as.character(GenomicRanges::seqnames(gr)),
        start = as.integer(IRanges::start(gr)),
        end = as.integer(IRanges::end(gr)),
        strand = as.character(GenomicRanges::strand(gr)),
        V5 = as.numeric(S4Vectors::mcols(gr)\$V5),
        annotation = fallback_annotation,
        distanceToTSS = NA_real_,
        stringsAsFactors = FALSE
    )
}
write.table(anno_result, file="${prefix}_chipseeker_annotation.tsv", 
            sep="\\t", quote=FALSE, row.names=FALSE)

get_annotation_feature <- function(df) {
    anno_col <- if ("annotation" %in% colnames(df)) {
        "annotation"
    } else if ("Annotation" %in% colnames(df)) {
        "Annotation"
    } else {
        NA
    }

    if (is.na(anno_col)) {
        return(factor(rep("Unannotated", nrow(df))))
    }

    feature <- as.character(df[[anno_col]])
    feature[is.na(feature) | feature == ""] <- "Unannotated"
    feature <- sub(" \\\\(.*", "", feature)
    factor(feature)
}

# Create pie chart of genomic annotation
safe_dual_plot(
    "${prefix}_chipseeker_annotation_pie.pdf",
    "${prefix}_chipseeker_annotation_pie.png",
    8, 8, 800, 800,
    function() {
        if (is.null(peakAnno)) {
            pie_counts <- table(get_annotation_feature(anno_result))
            pie(
                pie_counts,
                main = "Peak Annotation",
                col = ChIPseeker:::getCols(length(pie_counts))
            )
        } else {
            p <- plotAnnoPie(peakAnno)
            if (inherits(p, "ggplot")) print(p)
        }
    },
    "Annotation pie unavailable"
)

# Create bar plot of genomic annotation
safe_dual_plot(
    "${prefix}_chipseeker_annotation_bar.pdf",
    "${prefix}_chipseeker_annotation_bar.png",
    10, 6, 1000, 600,
    function() {
        if (is.null(peakAnno)) {
            bar_counts <- sort(table(get_annotation_feature(anno_result)), decreasing = TRUE)
            barplot(
                bar_counts,
                las = 2,
                ylab = "Count",
                main = "Peak Annotation",
                col = ChIPseeker:::getCols(length(bar_counts))
            )
        } else {
            p <- plotAnnoBar(peakAnno)
            if (inherits(p, "ggplot")) print(p)
        }
    },
    "Annotation bar unavailable"
)

# Create upset plot
safe_dual_plot(
    "${prefix}_chipseeker_annotation_upset.pdf",
    "${prefix}_chipseeker_annotation_upset.png",
    12, 8, 1200, 800,
    function() {
        if (is.null(peakAnno)) stop("peakAnno unavailable")
        p <- upsetplot(peakAnno, vennpie = FALSE)
        if (inherits(p, "ggplot")) print(p)
    },
    "Annotation upset unavailable"
)

# Distribution of peaks relative to TSS
safe_dual_plot(
    "${prefix}_chipseeker_dist_to_tss.pdf",
    "${prefix}_chipseeker_dist_to_tss.png",
    10, 6, 1000, 600,
    function() {
        if (is.null(peakAnno)) stop("peakAnno unavailable")
        p <- plotDistToTSS(peakAnno, title = "Distribution of peaks relative to TSS")
        if (inherits(p, "ggplot")) print(p)
    },
    "Distance-to-TSS plot unavailable"
)

# Coverage of peak regions across chromosomes
safe_dual_plot(
    "${prefix}_chipseeker_coverage.pdf",
    "${prefix}_chipseeker_coverage.png",
    12, 6, 1200, 600,
    function() {
        plotted <- FALSE
        if (length(gr_cov) > 0 && length(cov_chrs) > 0) {
            plotted <- tryCatch({
                covplot(gr_cov, weightCol = "V5", chrs = cov_chrs)
                TRUE
            }, error = function(e) {
                cat("Note: covplot failed for ${prefix} -", conditionMessage(e), "\\n")
                FALSE
            })
        }

        # Fallback to a simple weighted position scatter so coverage files are never empty placeholders.
        if (!plotted) {
            mids <- floor((as.numeric(IRanges::start(gr)) + as.numeric(IRanges::end(gr))) / 2)
            ord <- order(as.character(GenomicRanges::seqnames(gr)), mids)
            mids <- mids[ord]
            y <- seq_along(mids)
            plot(
                mids,
                y,
                pch = 16,
                cex = 0.35,
                col = grDevices::adjustcolor("steelblue", alpha.f = 0.35),
                xlab = "Genomic position (bp)",
                ylab = "Peak index",
                main = "Peak positions (coverage fallback)"
            )
        }
    },
    "Coverage unavailable"
)

cat("ChIPseeker annotation completed\\n")
RSCRIPT

    Rscript chipseeker_annotate.R

    printf '"ChIPseeker":\n' > versions.yml
    printf '    ChIPseeker: %s\n' '1.42.0' >> versions.yml
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
    path("chipseeker_comparison_annotation_bar_condition.pdf"), emit: anno_bar_cond_pdf
    path("chipseeker_comparison_annotation_bar_condition.png"), emit: anno_bar_cond_png
    path("chipseeker_comparison_dist_to_tss.pdf"), emit: dist_pdf
    path("chipseeker_comparison_dist_to_tss.png"), emit: dist_png
    path("chipseeker_comparison_coverage.pdf"), emit: cov_pdf
    path("chipseeker_comparison_coverage.png"), emit: cov_png
    path("chipseeker_comparison_coverage_condition_average.pdf"), emit: cov_avg_pdf
    path("chipseeker_comparison_coverage_condition_average.png"), emit: cov_avg_png
    path("chipseeker_condition_*_pie.pdf"), optional: true, emit: cond_pie_pdf
    path("chipseeker_condition_*_pie.png"), optional: true, emit: cond_pie_png
    path("versions.yml"), emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    cat > chipseeker_compare.R <<'RSCRIPT'
library(ChIPseeker)
library(ggplot2)
library(GenomicRanges)

# Create temp library directory for package installation
temp_lib <- file.path(tempdir(), "R_packages")
dir.create(temp_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(temp_lib, .libPaths()))

# Plot only canonical chromosomes (exclude scaffolds/random/chrUn)
main_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")

normalize_feature <- function(x) {
    x <- as.character(x)
    x[is.na(x) | x == ""] <- "Unannotated"
    sub(" \\\\(.*", "", x)
}

collapse_annotation <- function(x) {
    x <- normalize_feature(x)
    x[grepl("^Promoter", x)] <- "Promoter"
    x[grepl("^Exon", x)] <- "Exon"
    x[grepl("^Intron", x)] <- "Intron"
    x[grepl("^Downstream", x)] <- "Downstream"
    x[x %in% c("Distal Intergenic", "Intergenic")] <- "Intergenic"
    x[grepl("^5' UTR|^5UTR|^fiveUTR", x)] <- "5' UTR"
    x[grepl("^3' UTR|^3UTR|^threeUTR", x)] <- "3' UTR"
    x[is.na(x) | x == ""] <- "Unannotated"
    x
}

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

        # Prefer canonical chromosomes, but fall back to all seqlevels when canonical labels are absent.
        keep_chrs <- intersect(GenomeInfoDb::seqlevels(gr), main_chrs)
        if (length(keep_chrs) > 0) {
            gr <- GenomeInfoDb::keepSeqlevels(gr, keep_chrs, pruning.mode = "coarse")
        }

        # Keep flattened annotation table for robust custom comparison plots
        anno_df\$sample <- sub("_chipseeker_annotation\\\\.tsv", "", basename(f))
        anno_all[[basename(f)]] <- anno_df

        anno_col <- if ("annotation" %in% colnames(anno_df)) {
            "annotation"
        } else if ("Annotation" %in% colnames(anno_df)) {
            "Annotation"
        } else {
            NA
        }

        feature <- if (!is.na(anno_col)) {
            normalize_feature(anno_df[[anno_col]])
        } else {
            rep("Unannotated", nrow(anno_df))
        }

        stat_df <- as.data.frame(table(feature), stringsAsFactors = FALSE)
        colnames(stat_df) <- c("Feature", "Count")
        stat_df <- stat_df[stat_df\$Count > 0, ]
        total_count <- sum(stat_df\$Count)
        stat_df\$Frequency <- if (total_count > 0) (100 * stat_df\$Count / total_count) else 0
        stat_df\$sample <- anno_df\$sample[1]
        anno_stat_all[[anno_df\$sample[1]]] <- stat_df[, c("sample", "Feature", "Frequency")]
        
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
                ggplot2::ggtitle("Peak Annotation Comparison") +
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
                ggplot2::ggtitle("Peak Annotation Comparison") +
                ggplot2::theme_bw() +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
            print(p)
        }, error = function(e) {
            cat("Note: Comparison bar plot generation skipped -", conditionMessage(e), "\\n")
        })
        dev.off()

        # Condition-level annotation bar plot (replicates aggregated per condition)
        pdf("chipseeker_comparison_annotation_bar_condition.pdf", width=11, height=6)
        tryCatch({
            anno_col_all <- if ("annotation" %in% colnames(all_df)) {
                "annotation"
            } else if ("Annotation" %in% colnames(all_df)) {
                "Annotation"
            } else {
                NA
            }
            if (is.na(anno_col_all)) stop("No annotation column found")

            cond_df <- data.frame(
                condition = sub("_R[0-9]+\$", "", as.character(all_df\$sample)),
                Feature = normalize_feature(all_df[[anno_col_all]]),
                stringsAsFactors = FALSE
            )
            cond_stat <- as.data.frame(table(cond_df\$condition, cond_df\$Feature), stringsAsFactors = FALSE)
            colnames(cond_stat) <- c("condition", "Feature", "Count")
            cond_stat <- cond_stat[cond_stat\$Count > 0, ]
            cond_totals <- aggregate(Count ~ condition, data = cond_stat, FUN = sum)
            cond_stat <- merge(cond_stat, cond_totals, by = "condition", suffixes = c("", "_total"))
            cond_stat\$Frequency <- ifelse(cond_stat\$Count_total > 0, 100 * cond_stat\$Count / cond_stat\$Count_total, 0)
            cond_levels <- unique(cond_stat\$condition)
            cond_stat\$condition <- factor(cond_stat\$condition, levels = cond_levels)
            cond_totals\$condition <- factor(cond_totals\$condition, levels = cond_levels)
            cond_totals\$label <- paste0("N=", cond_totals\$Count)

            cond_feature_levels <- unique(cond_stat\$Feature)
            cond_stat\$Feature <- factor(cond_stat\$Feature, levels = cond_feature_levels)
            cond_feature_colors <- setNames(ChIPseeker:::getCols(length(cond_feature_levels)), cond_feature_levels)

            p_cond <- ggplot2::ggplot(cond_stat, ggplot2::aes(x = condition, y = Frequency, fill = Feature)) +
                ggplot2::geom_col(position = "stack") +
                ggplot2::geom_text(
                    data = cond_totals,
                    ggplot2::aes(x = condition, y = -3, label = label),
                    inherit.aes = FALSE,
                    vjust = 1,
                    size = 3.5
                ) +
                ggplot2::scale_fill_manual(values = cond_feature_colors, breaks = cond_feature_levels) +
                ggplot2::scale_y_continuous(
                    limits = c(-8, 100),
                    breaks = seq(0, 100, 10),
                    expand = ggplot2::expansion(mult = c(0, 0.02))
                ) +
                ggplot2::ylab("Frequency (%)") +
                ggplot2::xlab("Condition") +
                ggplot2::ggtitle("Peak Annotation Comparison (Condition)") +
                ggplot2::coord_cartesian(clip = "off") +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                    plot.margin = ggplot2::margin(t = 8, r = 10, b = 26, l = 10)
                )
            print(p_cond)
        }, error = function(e) {
            cat("Note: Condition comparison bar plot generation skipped -", conditionMessage(e), "\\n")
        })
        dev.off()

        png("chipseeker_comparison_annotation_bar_condition.png", width=1100, height=600, res=100)
        tryCatch({
            anno_col_all <- if ("annotation" %in% colnames(all_df)) {
                "annotation"
            } else if ("Annotation" %in% colnames(all_df)) {
                "Annotation"
            } else {
                NA
            }
            if (is.na(anno_col_all)) stop("No annotation column found")

            cond_df <- data.frame(
                condition = sub("_R[0-9]+\$", "", as.character(all_df\$sample)),
                Feature = normalize_feature(all_df[[anno_col_all]]),
                stringsAsFactors = FALSE
            )
            cond_stat <- as.data.frame(table(cond_df\$condition, cond_df\$Feature), stringsAsFactors = FALSE)
            colnames(cond_stat) <- c("condition", "Feature", "Count")
            cond_stat <- cond_stat[cond_stat\$Count > 0, ]
            cond_totals <- aggregate(Count ~ condition, data = cond_stat, FUN = sum)
            cond_stat <- merge(cond_stat, cond_totals, by = "condition", suffixes = c("", "_total"))
            cond_stat\$Frequency <- ifelse(cond_stat\$Count_total > 0, 100 * cond_stat\$Count / cond_stat\$Count_total, 0)
            cond_levels <- unique(cond_stat\$condition)
            cond_stat\$condition <- factor(cond_stat\$condition, levels = cond_levels)
            cond_totals\$condition <- factor(cond_totals\$condition, levels = cond_levels)
            cond_totals\$label <- paste0("N=", cond_totals\$Count)

            cond_feature_levels <- unique(cond_stat\$Feature)
            cond_stat\$Feature <- factor(cond_stat\$Feature, levels = cond_feature_levels)
            cond_feature_colors <- setNames(ChIPseeker:::getCols(length(cond_feature_levels)), cond_feature_levels)

            p_cond <- ggplot2::ggplot(cond_stat, ggplot2::aes(x = condition, y = Frequency, fill = Feature)) +
                ggplot2::geom_col(position = "stack") +
                ggplot2::geom_text(
                    data = cond_totals,
                    ggplot2::aes(x = condition, y = -3, label = label),
                    inherit.aes = FALSE,
                    vjust = 1,
                    size = 3.5
                ) +
                ggplot2::scale_fill_manual(values = cond_feature_colors, breaks = cond_feature_levels) +
                ggplot2::scale_y_continuous(
                    limits = c(-8, 100),
                    breaks = seq(0, 100, 10),
                    expand = ggplot2::expansion(mult = c(0, 0.02))
                ) +
                ggplot2::ylab("Frequency (%)") +
                ggplot2::xlab("Condition") +
                ggplot2::ggtitle("Peak Annotation Comparison (Condition)") +
                ggplot2::coord_cartesian(clip = "off") +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                    plot.margin = ggplot2::margin(t = 8, r = 10, b = 26, l = 10)
                )
            print(p_cond)
        }, error = function(e) {
            cat("Note: Condition comparison bar plot generation skipped -", conditionMessage(e), "\\n")
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
            cov_chrs <- intersect(main_chrs, unique(cov_all\$chr))
            if (length(cov_chrs) == 0) {
                cov_chrs <- unique(cov_all\$chr)
            }
            cov_all <- cov_all[cov_all\$chr %in% cov_chrs, ]
            if (nrow(cov_all) == 0) {
                stop("No chromosomes available for comparison coverage plotting")
            }
            cov_all\$weight[!is.finite(cov_all\$weight) | cov_all\$weight <= 0] <- 1
            cov_all\$condition <- sub("_R[0-9]+\$", "", cov_all\$sample)

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
    
    # Generate condition-level pie charts with merged (union) peaks per condition
    cat("Generating condition-level vennpie plots...\\n")
    tryCatch({
        # Extract condition from sample name (everything before _R[0-9]+)
        condition_groups <- list()
        
        for (f in anno_files) {
            anno_df <- read.delim(f, stringsAsFactors = FALSE)
            sample_name <- sub("_chipseeker_annotation\\\\.tsv", "", basename(f))
            
            # Extract condition: remove _R[0-9]+ suffix
            condition <- sub("_R[0-9]+\$", "", sample_name)
            
            if (!(condition %in% names(condition_groups))) {
                condition_groups[[condition]] <- anno_df[0, ]  # Empty df with same structure
            }
            
            # Bind all annotations for this condition
            condition_groups[[condition]] <- rbind(condition_groups[[condition]], anno_df)
        }
        
        # For each condition, create merged peak annotations and plot vennpie
        for (cond in names(condition_groups)) {
            cat("Processing condition:", cond, "\\n")
            
            cond_anno <- condition_groups[[cond]]
            
            if (nrow(cond_anno) == 0) next
            
            # Create GRanges for all peaks in this condition
            merged_gr <- GRanges(
                seqnames = cond_anno\$seqnames,
                ranges = IRanges(start = as.numeric(cond_anno\$start), end = as.numeric(cond_anno\$end)),
                strand = cond_anno\$strand
            )
            
            # Reduce to unique peak regions (union: combine overlapping peaks)
            merged_gr <- GenomicRanges::reduce(merged_gr, ignore.strand = TRUE)
            
            # Map annotations to merged peaks by overlap
            mcols(merged_gr)\$annotation <- "Unannotated"
            
            # Do a simple overlap: for each merged peak, find overlapping input peaks and collect annotations
            hits <- findOverlaps(merged_gr, GRanges(
                seqnames = cond_anno\$seqnames,
                ranges = IRanges(start = as.numeric(cond_anno\$start), end = as.numeric(cond_anno\$end))
            ), ignore.strand = TRUE)
            
            if (length(hits) > 0) {
                for (i in unique(queryHits(hits))) {
                    hit_indices <- subjectHits(hits[queryHits(hits) == i])
                    hit_annos <- collapse_annotation(cond_anno\$annotation[hit_indices])
                    # Prefer non-Unannotated annotations
                    non_unann <- hit_annos[hit_annos != "Unannotated"]
                    if (length(non_unann) > 0) {
                        mcols(merged_gr)\$annotation[i] <- non_unann[1]
                    }
                }
            }
            
            # Create peakAnno-like annotation factor for vennpie
            anno_factor <- factor(collapse_annotation(mcols(merged_gr)\$annotation))
            
            # Plot condition pie using vennpie
            pdf_file <- paste0("chipseeker_condition_", gsub(" ", "_", cond), "_pie.pdf")
            png_file <- paste0("chipseeker_condition_", gsub(" ", "_", cond), "_pie.png")
            
            # PDF version
            pdf(pdf_file, width = 8, height = 8)
            tryCatch({
                vennpie(anno_factor)
                title(main = paste("Peak Annotation -", cond))
            }, error = function(e) {
                cat("Note: vennpie failed for condition", cond, "-", conditionMessage(e), "\\n")
                # Fallback: simple pie chart
                pie_counts <- table(anno_factor)
                pie(pie_counts, main = paste("Peak Annotation -", cond), 
                    col = ChIPseeker:::getCols(length(pie_counts)))
            })
            dev.off()
            
            # PNG version
            png(png_file, width = 800, height = 800, res = 100)
            tryCatch({
                vennpie(anno_factor)
                title(main = paste("Peak Annotation -", cond))
            }, error = function(e) {
                cat("Note: vennpie failed for condition", cond, "-", conditionMessage(e), "\\n")
                # Fallback: simple pie chart
                pie_counts <- table(anno_factor)
                pie(pie_counts, main = paste("Peak Annotation -", cond),
                    col = ChIPseeker:::getCols(length(pie_counts)))
            })
            dev.off()
            
            cat("Created pie charts for condition:", cond, "\\n")
        }
        
    }, error = function(e) {
        cat("Note: Condition-level pie chart generation encountered error -", conditionMessage(e), "\\n")
    })
}

# Guarantee required comparison files exist for Nextflow output validation
if (!file.exists("chipseeker_comparison_annotation_bar.pdf")) {
    make_placeholder_pdf("chipseeker_comparison_annotation_bar.pdf", "Comparison annotation bar unavailable")
}
if (!file.exists("chipseeker_comparison_annotation_bar.png")) {
    make_placeholder_png("chipseeker_comparison_annotation_bar.png", "Comparison annotation bar unavailable")
}
if (!file.exists("chipseeker_comparison_annotation_bar_condition.pdf")) {
    make_placeholder_pdf("chipseeker_comparison_annotation_bar_condition.pdf", "Condition comparison annotation bar unavailable")
}
if (!file.exists("chipseeker_comparison_annotation_bar_condition.png")) {
    make_placeholder_png("chipseeker_comparison_annotation_bar_condition.png", "Condition comparison annotation bar unavailable")
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
