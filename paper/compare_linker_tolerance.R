#!/usr/bin/env Rscript
## Compares runs of one sample that differ in cb_umi_max_errors: cells, median
## UMIs per cell, and pseudobulk agreement with the first run, per aligner.
##
## Rscript paper/compare_linker_tolerance.R \
##     --runs 0=output/sendoel2024 1=output/sendoel2024_linker1 2=output/sendoel2024_linker2 \
##     --sample sample_16_wta_p60 --aligners starsolo,kallisto,alevin \
##     --out_prefix output/sendoel2024/paper/sample_16_wta_p60/linker_tolerance

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(ggplot2)
})

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(script_dir)) script_dir <- "paper"
source(file.path(script_dir, "..", "workflow", "src", "agreement_metrics.R"))

run_summary <- function(counts) {
    data.table(n_cells = ncol(counts), median_umis = median(Matrix::colSums(counts)))
}

read_counts <- function(working_dir, aligner, sample) {
    path <- file.path(working_dir, aligner, sample, sprintf("%s_%s_sce.rds", sample, aligner))
    SingleCellExperiment::counts(readRDS(path))
}

main <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--runs", required = TRUE, nargs = "+",
                        help = "max_errors=working_dir pairs; the first is the reference")
    parser$add_argument("--sample", required = TRUE)
    parser$add_argument("--aligners", required = TRUE, help = "comma-separated")
    parser$add_argument("--out_prefix", required = TRUE)
    args <- parser$parse_args()

    runs <- strsplit(args$runs, "=", fixed = TRUE)
    aligners <- strsplit(args$aligners, ",", fixed = TRUE)[[1]]

    comparison <- rbindlist(lapply(aligners, function(aligner) {
        counts <- lapply(runs, function(run) read_counts(run[2], aligner, args$sample))
        pseudobulk <- lapply(counts, Matrix::rowSums)
        rbindlist(lapply(seq_along(runs), function(i) {
            cbind(aligner = aligner, max_errors = runs[[i]][1], run_summary(counts[[i]]),
                  pseudobulk_agreement(pseudobulk[[1]], pseudobulk[[i]]))
        }))
    }))

    dir.create(dirname(args$out_prefix), recursive = TRUE, showWarnings = FALSE)
    fwrite(comparison, paste0(args$out_prefix, ".csv"))
    long <- melt(comparison, id.vars = c("aligner", "max_errors"),
                 measure.vars = c("n_cells", "median_umis", "scaled_mard"))
    plot <- ggplot(long, aes(max_errors, value, colour = aligner, group = aligner)) +
        geom_line() +
        geom_point() +
        facet_wrap(~variable, scales = "free_y") +
        theme_bw() +
        labs(x = "cb_umi_max_errors", y = NULL, colour = "aligner")
    ggsave(paste0(args$out_prefix, ".pdf"), plot, width = 8, height = 2.8)
    print(comparison)
}

if (sys.nframe() == 0) main()
