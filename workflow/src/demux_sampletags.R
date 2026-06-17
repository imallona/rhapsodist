#!/usr/bin/env R
##
## Shared sampletag demultiplexing. Reads the per-read count table (cb, umi,
## sampletag, extra), builds a cell by tag matrix and runs the BD algorithm:
## high-quality singlets (>75% of reads from one tag), a linear noise estimate,
## then per-tag minimum-count thresholds for the rest. allowed_tags restricts the
## candidate tags first (covers both the starsolo and search paths). Sourced by
## generate_sampletag_report.Rmd and the splitting rule so the assignment is shared.
##
## GPLv3

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
})

demux_sampletags <- function(counts_path, allowed_tags = NULL) {
    counts_raw <- fread(
        cmd = paste("zcat", shQuote(counts_path)),
        sep = "\t", header = FALSE,
        col.names = c("cb", "umi", "sampletag", "extra")
    )

    agg <- counts_raw[, .N, by = .(cb, sampletag)]
    cb_f  <- factor(agg$cb)
    tag_f <- factor(agg$sampletag)
    mat <- as.matrix(sparseMatrix(
        i = as.integer(cb_f),
        j = as.integer(tag_f),
        x = agg$N,
        dims = c(nlevels(cb_f), nlevels(tag_f)),
        dimnames = list(levels(cb_f), levels(tag_f))
    ))

    if (!is.null(allowed_tags)) {
        keep <- intersect(colnames(mat), allowed_tags)
        mat <- mat[, keep, drop = FALSE]
    }

    # step 1: high-quality singlets
    total <- rowSums(mat)
    props <- mat / pmax(total, 1L)
    dominant_idx <- max.col(props, ties.method = "first")
    dominant_frac <- props[cbind(seq_len(nrow(props)), dominant_idx)]
    dominant_tag <- colnames(mat)[dominant_idx]
    is_highqual <- dominant_frac >= 0.75 & total > 0

    # step 2: per-tag thresholds (min count among high-quality singlets for each tag)
    tag_thresholds <- vapply(seq_along(colnames(mat)), function(j) {
        hq <- is_highqual & dominant_idx == j
        if (!any(hq)) return(NA_real_)
        min(mat[hq, j])
    }, numeric(1))
    names(tag_thresholds) <- colnames(mat)

    # step 3: noise model using high-quality singlets only
    hq_counts <- mat[cbind(seq_len(nrow(mat)), dominant_idx)]
    noise_hq <- ifelse(is_highqual, total - hq_counts, NA_real_)

    # per-tag noise fraction: fraction of background attributable to each tag
    noise_mat <- mat
    noise_mat[cbind(which(is_highqual), dominant_idx[is_highqual])] <- 0L
    noise_mat[!is_highqual, ] <- NA_real_
    tag_noise_frac <- colSums(noise_mat, na.rm = TRUE)
    total_noise <- sum(tag_noise_frac)
    if (total_noise > 0) {
        tag_noise_frac <- tag_noise_frac / total_noise
    } else {
        tag_noise_frac[] <- 1 / ncol(mat)
    }

    # analytic ols slope (no intercept): noise ~ 0 + total
    hq_total <- total[is_highqual]
    hq_noise <- noise_hq[is_highqual]
    slope <- if (length(hq_total) >= 2 && var(hq_total) > 0) {
        max(sum(hq_total * hq_noise) / sum(hq_total^2), 0)
    } else {
        0
    }

    # step 4: noise-corrected counts
    # corrected[i, j] = mat[i, j] - slope * total[i] * tag_noise_frac[j]
    predicted_noise <- slope * total
    corrected <- mat - outer(predicted_noise, tag_noise_frac)

    # step 5: call non-highqual cells against per-tag thresholds
    valid <- !is.na(tag_thresholds)
    called_mat <- matrix(FALSE, nrow(mat), ncol(mat), dimnames = dimnames(mat))
    if (any(valid)) {
        called_mat[, valid] <- sweep(
            corrected[, valid, drop = FALSE], 2, tag_thresholds[valid], ">"
        )
    }
    n_called <- rowSums(called_mat)

    # final classification
    status <- rep("undetermined", nrow(mat))
    called_tag <- rep(NA_character_, nrow(mat))

    status[is_highqual] <- "highqual"
    called_tag[is_highqual] <- dominant_tag[is_highqual]

    is_called <- !is_highqual & n_called == 1
    status[is_called] <- "called"
    if (any(is_called)) {
        called_tag[is_called] <- colnames(called_mat)[
            max.col(called_mat[is_called, , drop = FALSE], ties.method = "first")
        ]
    }
    status[!is_highqual & n_called > 1] <- "multiplet"

    demux <- data.frame(
        cb = rownames(mat),
        called_tag = called_tag,
        dominant_frac = round(dominant_frac, 3),
        total_reads = total,
        noise_reads = noise_hq,
        predicted_noise = round(predicted_noise, 1),
        status = status,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    attr(demux, "slope") <- slope
    demux
}

## CLI: only runs when invoked as a script with --counts (the Rmd sources this file
## to obtain demux_sampletags() and must not trigger the CLI).
.cli_args <- commandArgs(trailingOnly = TRUE)
if (length(.cli_args) > 0 && "--counts" %in% .cli_args) {
    suppressPackageStartupMessages(library(argparse))
    parser <- ArgumentParser(description = "Demultiplex sampletag counts into a per-cell assignment table.")
    parser$add_argument("--counts", required = TRUE, help = "per-read sampletag counts tsv.gz")
    parser$add_argument("--out", required = TRUE, help = "output demux tsv.gz")
    parser$add_argument("--rds", default = "", help = "optional demux rds for debugging")
    parser$add_argument("--allowed_tags", default = "",
                        help = "optional comma-separated tag names to restrict candidates")
    args <- parser$parse_args()

    allowed <- if (nzchar(args$allowed_tags)) strsplit(args$allowed_tags, ",")[[1]] else NULL
    demux <- demux_sampletags(args$counts, allowed_tags = allowed)

    cat(sprintf("noise model slope: %.4f\n", attr(demux, "slope")))
    print(table(demux$status))
    print(table(demux$called_tag, useNA = "ifany"))

    fwrite(demux, args$out, sep = "\t", compress = "gzip")
    if (nzchar(args$rds)) saveRDS(demux, args$rds)
}
