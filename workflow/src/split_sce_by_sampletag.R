#!/usr/bin/env R
##
## Split an aligner SingleCellExperiment into one object per sampletag, using the
## demux_sampletags.R assignment. For each tag in use, writes the singlet cells
## (high-quality and called) as an HDF5-backed SCE via saveHDF5SummarizedExperiment;
## no top-level *_sce.rds, so the descriptive report's sce.rds glob skips them.
## --tags gives the configured tags (empty objects included for tags with no cells);
## otherwise the observed tags are used.
##
## GPLv3

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(HDF5Array)
    library(data.table)
    library(argparse)
})

parser <- ArgumentParser(description = "Split an aligner SCE into per-sampletag HDF5 SCEs.")
parser$add_argument("--sce", required = TRUE, help = "aligner SCE rds (HDF5-backed)")
parser$add_argument("--demux", required = TRUE, help = "demux assignment tsv.gz")
parser$add_argument("--out_base", required = TRUE, help = "base directory for by_sampletag output")
parser$add_argument("--tags", default = "", help = "optional comma-separated tag names to produce")
parser$add_argument("--labels", default = "", help = "optional comma-separated output labels, parallel to --tags")
parser$add_argument("--summary_rds", default = "", help = "optional summary rds for debugging")
parser$add_argument("--done", default = "", help = "optional sentinel file to touch when finished")
args <- parser$parse_args()

sce <- readRDS(args$sce)
demux <- as.data.frame(fread(cmd = paste("zcat", shQuote(args$demux)), sep = "\t", header = TRUE))

singlets <- demux[demux$status %in% c("highqual", "called") & !is.na(demux$called_tag), ]
overlap <- length(intersect(singlets$cb, colnames(sce)))
cat(sprintf("singlet cells: %d; present in SCE: %d (%.1f%%)\n",
            nrow(singlets), overlap,
            100 * overlap / max(nrow(singlets), 1)))

if (nzchar(args$tags)) {
    tags <- strsplit(args$tags, ",")[[1]]
    labels <- if (nzchar(args$labels)) strsplit(args$labels, ",")[[1]] else tags
    stopifnot(length(tags) == length(labels))
} else {
    tags <- sort(unique(singlets$called_tag))
    labels <- tags
}

dir.create(args$out_base, recursive = TRUE, showWarnings = FALSE)
summary <- data.frame(tag = character(), label = character(), n_cells = integer(),
                      stringsAsFactors = FALSE)

for (i in seq_along(tags)) {
    barcodes <- singlets$cb[singlets$called_tag == tags[i]]
    cells <- intersect(barcodes, colnames(sce))
    sub <- sce[, cells, drop = FALSE]
    out_dir <- file.path(args$out_base, labels[i])
    saveHDF5SummarizedExperiment(sub, dir = out_dir, replace = TRUE)
    cat(sprintf("tag %s -> %s: %d cells\n", tags[i], labels[i], length(cells)))
    summary <- rbind(summary, data.frame(tag = tags[i], label = labels[i],
                                          n_cells = length(cells), stringsAsFactors = FALSE))
}

if (nzchar(args$summary_rds)) saveRDS(summary, args$summary_rds)
if (nzchar(args$done)) file.create(args$done)
