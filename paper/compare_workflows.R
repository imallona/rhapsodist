#!/usr/bin/env Rscript
## Compares other workflows with the rhapsodist aligners on one sample: cells,
## pseudobulk agreement with each aligner, time and memory. A workflow without
## an SCE is listed as failed.
##
## Rscript paper/compare_workflows.R \
##     --working_dir output/sendoel2024 --sample sample_16_wta_p60 --bead_version v1 \
##     --aligners starsolo,kallisto,alevin,sbg --workflows universc,zumis,openpipelines \
##     --out_prefix output/sendoel2024/paper/sample_16_wta_p60/workflow_comparison

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(ggplot2)
})

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(script_dir)) script_dir <- "paper"
source(file.path(script_dir, "..", "workflow", "src", "agreement_metrics.R"))

sce_path <- function(working_dir, tool, sample) {
    file.path(working_dir, tool, sample, sprintf("%s_%s_sce.rds", sample, tool))
}

## Seconds and peak memory of a Snakemake benchmark file; NA when it is absent.
read_benchmark <- function(path) {
    if (!file.exists(path)) {
        return(data.table(seconds = NA_real_, max_rss_mb = NA_real_))
    }
    benchmark <- fread(path)
    data.table(seconds = sum(benchmark$s), max_rss_mb = max(benchmark$max_rss))
}

## Barcodes and pseudobulk of a count matrix, all that the comparison needs.
tool_summary <- function(counts) {
    list(cells = colnames(counts), pseudobulk = Matrix::rowSums(counts))
}

## One row per workflow and aligner. summary is NULL for a failed workflow.
workflow_rows <- function(workflow, summary, aligner_summaries, benchmark) {
    rows <- rbindlist(lapply(names(aligner_summaries), function(aligner) {
        if (is.null(summary)) {
            return(data.table(aligner = aligner, shared_cells = NA_integer_,
                              spearman = NA_real_, scaled_mard = NA_real_))
        }
        reference <- aligner_summaries[[aligner]]
        cbind(data.table(aligner = aligner,
                         shared_cells = length(intersect(summary$cells, reference$cells))),
              pseudobulk_agreement(reference$pseudobulk, summary$pseudobulk))
    }))
    cbind(data.table(workflow = workflow,
                     status = if (is.null(summary)) "failed" else "ok",
                     n_cells = if (is.null(summary)) NA_integer_ else length(summary$cells)),
          rows, benchmark)
}

main <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--working_dir", required = TRUE)
    parser$add_argument("--sample", required = TRUE)
    parser$add_argument("--bead_version", required = TRUE, help = "written to the table")
    parser$add_argument("--aligners", required = TRUE, help = "comma-separated")
    parser$add_argument("--workflows", required = TRUE, help = "comma-separated")
    parser$add_argument("--out_prefix", required = TRUE)
    args <- parser$parse_args()

    ## the HDF5-backed counts are read into memory once per tool
    read_summary <- function(tool) {
        sce <- readRDS(sce_path(args$working_dir, tool, args$sample))
        tool_summary(as(SingleCellExperiment::counts(sce), "CsparseMatrix"))
    }
    aligners <- strsplit(args$aligners, ",", fixed = TRUE)[[1]]
    aligner_summaries <- lapply(setNames(aligners, aligners), read_summary)

    comparison <- rbindlist(lapply(strsplit(args$workflows, ",", fixed = TRUE)[[1]], function(workflow) {
        has_sce <- file.exists(sce_path(args$working_dir, workflow, args$sample))
        benchmark <- read_benchmark(file.path(args$working_dir, "benchmarks",
                                              sprintf("%s_%s.txt", workflow, args$sample)))
        workflow_rows(workflow, if (has_sce) read_summary(workflow) else NULL,
                      aligner_summaries, benchmark)
    }))
    comparison <- cbind(sample = args$sample, bead_version = args$bead_version, comparison)

    dir.create(dirname(args$out_prefix), recursive = TRUE, showWarnings = FALSE)
    fwrite(comparison, paste0(args$out_prefix, ".csv"))

    cells <- rbind(unique(comparison[status == "ok", .(tool = workflow, n_cells)]),
                   data.table(tool = aligners,
                              n_cells = vapply(aligner_summaries, function(s) length(s$cells),
                                               integer(1))))
    cells_plot <- ggplot(cells, aes(tool, n_cells)) +
        geom_col() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        labs(x = NULL, y = "cells")
    ggsave(paste0(args$out_prefix, "_cells.pdf"), cells_plot, width = 4, height = 3)
    agreement_plot <- ggplot(comparison[status == "ok"], aes(aligner, workflow, fill = spearman)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%.3f", spearman)), size = 3) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        theme_bw() +
        labs(x = "rhapsodist aligner", y = NULL, fill = "Spearman")
    ggsave(paste0(args$out_prefix, "_agreement.pdf"), agreement_plot, width = 5, height = 2.5)
    print(comparison)
}

if (sys.nframe() == 0) main()
