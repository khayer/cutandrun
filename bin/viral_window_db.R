#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# viral_window_db.R
# Authors : Katharina (Nina) Hayer, co-created with Claude Opus 5 (Anthropic)
# Version : 1.4.0
# Date    : 2026-11-25
#
# Purpose : Differential coverage on a small viral contig, the csaw way:
#           fixed-width windows, counts per sample, a negative-binomial fit
#           whose dispersion is estimated FROM THE REPLICATES, then adjacent
#           significant windows merged and their p-values combined by Simes.
#
# Why not MACS2: on a 36 kb contig at ~24,000x, target and control both
#           saturate and a peak caller has no model of how much coverage varies
#           between replicates by chance. Measured on the C14 Ad5 data, two IgG
#           replicates produced 23 peaks while the real comparisons produced
#           15-21, and 14 of 18 "reproducible" target peaks were also called by
#           IgG vs IgG. The replicates carry the information the peak caller
#           never sees. See claude/weitzman-jack/ad5_peak_calling_limit.md.
#
# The same control-vs-control test is built in and always reported, because the
# failure above is invisible without it.
#
# Depends on edgeR and limma only. The csaw package is not required: the
# windowing is done upstream with bedtools and the merging is implemented here,
# so the statistics are visible in this file rather than inside a package.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
    library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
    i <- match(flag, args)
    if (is.na(i)) return(default)
    args[i + 1]
}
has <- function(flag) flag %in% args

counts_file  <- getarg("--counts")
samples_file <- getarg("--samples")
out_prefix   <- getarg("--out-prefix", "viral_db")
fdr_cut      <- as.numeric(getarg("--fdr", "0.05"))
merge_gap    <- as.numeric(getarg("--merge-gap", "200"))
min_count    <- as.numeric(getarg("--min-count", "10"))
norm_mode    <- getarg("--norm", "tmm")          # tmm | spikein | none
contrasts_f  <- getarg("--contrasts")            # optional TSV: treatment<TAB>control
                                                 # same format as --deseq2_contrasts

if (is.null(counts_file) || is.null(samples_file)) {
    stop("need --counts and --samples")
}

counts_tab <- read.delim(counts_file, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
samples    <- read.delim(samples_file, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)

coord_cols <- c("chrom", "start", "end")
if (!all(coord_cols %in% colnames(counts_tab))) {
    stop("counts file must begin with columns chrom, start, end")
}
coords  <- counts_tab[, coord_cols]
sample_cols <- setdiff(colnames(counts_tab), coord_cols)

if (!all(samples$sample %in% sample_cols)) {
    missing <- setdiff(samples$sample, sample_cols)
    stop("samples not present as count columns: ", paste(missing, collapse = ", "))
}
mat <- as.matrix(counts_tab[, samples$sample, drop = FALSE])
storage.mode(mat) <- "numeric"

message("windows: ", nrow(mat), "   samples: ", ncol(mat))

# -----------------------------------------------------------------------------
# One test: a set of "treatment" samples against a set of "reference" samples.
# Returns merged regions with a Simes-combined p-value, and the window table.
# -----------------------------------------------------------------------------
run_test <- function(mat, coords, grp, label) {
    keep_samples <- !is.na(grp)
    m   <- mat[, keep_samples, drop = FALSE]
    g   <- factor(grp[keep_samples], levels = c("ref", "trt"))

    if (length(unique(g)) < 2) {
        return(list(ok = FALSE, reason = "only one condition present"))
    }
    if (min(table(g)) < 2) {
        return(list(ok = FALSE,
                    reason = paste0("fewer than 2 replicates in a condition (",
                                    paste(names(table(g)), table(g), sep = "=",
                                          collapse = ", "), ")")))
    }

    y <- DGEList(counts = m, group = g)

    if (norm_mode == "spikein") {
        if (!"scale_factor" %in% colnames(samples)) {
            stop("--norm spikein needs a scale_factor column in the samples file")
        }
        sf <- samples$scale_factor[match(colnames(m), samples$sample)]
        if (any(is.na(sf)) || any(sf <= 0)) stop("bad scale_factor values")
        # edgeR divides by norm.factors, so a larger scale factor (a sample that
        # needs scaling UP) must become a smaller norm factor.
        nf <- 1 / sf
        y$samples$norm.factors <- nf / exp(mean(log(nf)))
    } else if (norm_mode == "tmm") {
        y <- calcNormFactors(y, method = "TMM")
    } else {
        y$samples$norm.factors <- rep(1, ncol(m))
    }

    # Keep windows with real coverage in at least as many samples as the
    # smallest condition, so a window carried by one library cannot drive a call.
    keep <- rowSums(m >= min_count) >= min(table(g))
    if (sum(keep) == 0) {
        return(list(ok = FALSE, reason = "no window passed the count filter"))
    }
    y <- y[keep, , keep.lib.sizes = FALSE]
    win <- coords[keep, , drop = FALSE]

    design <- model.matrix(~ g)
    y   <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    res <- glmQLFTest(fit, coef = 2)          # trt vs ref

    tab <- data.frame(win,
                      logFC  = res$table$logFC,
                      logCPM = res$table$logCPM,
                      PValue = res$table$PValue,
                      stringsAsFactors = FALSE)
    tab$FDR_window <- p.adjust(tab$PValue, method = "BH")

    # --- regions are runs of SIGNIFICANT windows ------------------------------
    # The windows tile the contig, so clustering every window by proximity
    # collapses the whole genome into one region - the first version of this
    # script did exactly that and reported "1 region" spanning 0-35000. FDR is
    # controlled at the window level by BH across all tested windows; a region
    # is simply a run of adjacent windows that passed, and its Simes p-value is
    # descriptive rather than a second controlled quantity.
    ord <- order(tab$chrom, tab$start)
    tab <- tab[ord, ]
    sig <- which(tab$FDR_window < fdr_cut)

    simes <- function(p) {
        p <- sort(p)
        min(1, min(length(p) * p / seq_along(p)))
    }

    if (length(sig) == 0) {
        regions <- data.frame(chrom = character(0), start = numeric(0),
                              end = numeric(0), n_windows = integer(0),
                              best_logFC = numeric(0), max_logCPM = numeric(0),
                              simes_p = numeric(0), min_FDR_window = numeric(0),
                              stringsAsFactors = FALSE)
    } else {
        st  <- tab[sig, ]
        brk <- c(TRUE, !(st$chrom[-1] == st$chrom[-nrow(st)] &
                         (st$start[-1] - st$end[-nrow(st)]) <= merge_gap))
        st$region <- cumsum(brk)
        regions <- do.call(rbind, lapply(split(seq_len(nrow(st)), st$region), function(i) {
            sub  <- st[i, ]
            best <- which.min(sub$PValue)
            data.frame(chrom          = sub$chrom[1],
                       start          = min(sub$start),
                       end            = max(sub$end),
                       n_windows      = nrow(sub),
                       best_logFC     = sub$logFC[best],
                       max_logCPM     = max(sub$logCPM),
                       simes_p        = simes(sub$PValue),
                       min_FDR_window = min(sub$FDR_window),
                       stringsAsFactors = FALSE)
        }))
        regions <- regions[order(regions$min_FDR_window, -abs(regions$best_logFC)), ]
    }

    list(ok = TRUE, windows = tab, regions = regions,
         n_sig_windows = length(sig), n_regions = nrow(regions),
         n_tested_windows = nrow(tab), label = label)
}

write_result <- function(r, prefix) {
    if (!isTRUE(r$ok)) return(invisible(NULL))
    write.table(r$windows, paste0(prefix, ".windows.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    write.table(r$regions, paste0(prefix, ".regions.tsv"), sep = "\t",
                quote = FALSE, row.names = FALSE)
    sig <- r$regions
    if (nrow(sig) > 0) {
        bed <- data.frame(sig$chrom, sig$start, sig$end,
                          sprintf("region_%d", seq_len(nrow(sig))),
                          round(-10 * log10(pmax(sig$min_FDR_window, 1e-300))), ".")
        write.table(bed, paste0(prefix, ".significant.bed"), sep = "\t",
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
        invisible(file.create(paste0(prefix, ".significant.bed")))
    }
}

# -----------------------------------------------------------------------------
# One contrast per target group: that group against the control group its
# samplesheet names. Pooling all targets against all controls would put the
# timepoint difference inside the contrast, which is not what is being asked.
# -----------------------------------------------------------------------------
is_ctl <- samples$is_control %in% c(1, "1", TRUE, "true", "TRUE")
if (!"control_group" %in% colnames(samples)) {
    stop("samples file needs a control_group column naming each target's control")
}
# Two kinds of question, both answered the same way:
#   no --contrasts : each target group against the control group its samplesheet
#                    names. "Is there any enrichment here at all?"
#   --contrasts    : the pairs listed in the file, e.g. ADNP_B02_25 vs
#                    ADNP_DMSO_25. "Does it change between conditions?" This is
#                    the question run_deseq2_peak_analysis asks, but over fixed
#                    windows rather than over consensus peaks - which matters on
#                    a viral contig, where the peaks themselves are not reliable.
if (!is.null(contrasts_f)) {
    ct <- read.delim(contrasts_f, header = FALSE, comment.char = "#",
                     stringsAsFactors = FALSE)
    ct <- ct[nzchar(trimws(ct[[1]])), , drop = FALSE]
    if (ncol(ct) < 2) stop("--contrasts needs two tab-separated columns")
    pairs <- data.frame(trt = trimws(ct[[1]]), ref = trimws(ct[[2]]),
                        stringsAsFactors = FALSE)
} else {
    tg <- unique(samples$group[!is_ctl])
    pairs <- data.frame(
        trt = tg,
        ref = vapply(tg, function(g) {
            cg <- unique(samples$control_group[samples$group == g])
            cg <- cg[!is.na(cg) & nzchar(cg)]
            if (length(cg) == 1) cg else NA_character_
        }, character(1)),
        stringsAsFactors = FALSE)
}

# Permutation null. With 3 replicates a side there is no way to split one
# control group in two and still estimate dispersion, so the null is built by
# relabelling: every balanced partition of the same samples that is not the
# true one. If an arbitrary labelling finds as much as the real one does, the
# design cannot tell occupancy from replicate variation.
balanced_partitions <- function(n_total, n_trt, truth, max_n = 50) {
    combos <- combn(n_total, n_trt, simplify = FALSE)
    seen <- list(); out <- list()
    comp_truth <- setdiff(seq_len(n_total), truth)
    for (cb in combos) {
        key <- paste(sort(cb), collapse = ",")
        ckey <- paste(sort(setdiff(seq_len(n_total), cb)), collapse = ",")
        if (key %in% seen || ckey %in% seen) next
        seen <- c(seen, key)
        # The complement of the true split is the SAME partition with the labels
        # flipped, and the test is two-sided, so it reproduces the real result
        # exactly. Excluding only the truth leaves it in the null.
        if (setequal(cb, truth) || setequal(cb, comp_truth)) next
        out <- c(out, list(cb))
    }
    if (length(out) > max_n) out <- out[sample(seq_along(out), max_n)]
    out
}

con <- file(paste0(out_prefix, ".summary.txt"), "w")
emit <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }
emit("viral_window_db.R 1.4.0   FDR < ", fdr_cut, "   norm = ", norm_mode,
     "   merge gap ", merge_gap, " bp")
# edgeR 4.4 changed the default quasi-likelihood method (legacy = FALSE), so
# results are not identical across that boundary. Record which one ran.
emit("edgeR ", as.character(packageVersion("edgeR")),
     "   R ", paste0(R.version$major, ".", R.version$minor))
emit("windows supplied: ", nrow(mat), "   samples: ", ncol(mat))
emit("")

any_run <- FALSE
for (k in seq_len(nrow(pairs))) {
    tg <- pairs$trt[k]
    cg <- pairs$ref[k]
    if (is.na(cg) || !nzchar(cg)) {
        emit(tg, ": SKIPPED - no control group named for it")
        next
    }
    idx_t <- which(samples$group == tg)
    if (length(idx_t) == 0) {
        emit(tg, ": SKIPPED - no samples found for group ", tg)
        next
    }
    idx_c <- which(samples$group == cg)
    if (length(idx_c) == 0) {
        emit(tg, ": SKIPPED - no samples found for control group ", cg)
        next
    }
    use <- c(idx_c, idx_t)
    grp <- rep(NA_character_, nrow(samples))
    grp[idx_c] <- "ref"; grp[idx_t] <- "trt"

    real <- run_test(mat, coords, grp, paste0(tg, " vs ", cg))
    if (!isTRUE(real$ok)) {
        emit(tg, " vs ", cg, ": NOT RUN - ", real$reason)
        next
    }
    any_run <- TRUE
    prefix <- paste0(out_prefix, ".", tg, "_vs_", cg)
    write_result(real, prefix)

    # --- permutation null on exactly these samples ---------------------------
    truth <- which(use %in% idx_t)
    parts <- balanced_partitions(length(use), length(idx_t), truth,
                                 max_n = as.numeric(getarg("--max-permutations", "50")))
    perm_counts <- integer(0)
    for (cb in parts) {
        g2 <- rep(NA_character_, nrow(samples))
        g2[use] <- "ref"; g2[use[cb]] <- "trt"
        r2 <- run_test(mat, coords, g2, "permuted")
        if (isTRUE(r2$ok)) perm_counts <- c(perm_counts, r2$n_sig_windows)
    }
    write.table(data.frame(permutation = seq_along(perm_counts),
                           significant_windows = perm_counts),
                paste0(prefix, ".permutation_null.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    emit(tg, " vs ", cg, ":")
    emit("  real          : ", real$n_sig_windows, " significant windows in ",
         real$n_regions, " regions (of ", real$n_tested_windows, " tested)")
    if (length(perm_counts) == 0) {
        emit("  permuted null : NOT RUN - no balanced relabelling was testable")
        emit("  WARNING: without a null the count above cannot be interpreted.")
    } else {
        emp_p <- (1 + sum(perm_counts >= real$n_sig_windows)) / (1 + length(perm_counts))
        emit("  permuted null : ", length(perm_counts), " relabellings, significant windows ",
             min(perm_counts), "-", max(perm_counts),
             " (median ", stats::median(perm_counts), ")")
        emit("  empirical p   : ", signif(emp_p, 3),
             "  (floor ", signif(1 / (1 + length(perm_counts)), 3),
             " with this many relabellings)")
        if (real$n_sig_windows == 0) {
            emit("  Nothing called, and the null agrees. No enrichment detectable here.")
        } else if (real$n_sig_windows <= max(perm_counts)) {
            emit("  WARNING: an arbitrary relabelling of these same samples reaches the")
            emit("  real result. This is replicate variation, not occupancy.")
        } else {
            emit("  The real labelling beats every relabelling tried.")
        }
    }
    emit("")
}
if (!any_run) emit("No contrast could be tested.")
close(con)
