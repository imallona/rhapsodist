#!/usr/bin/env R
##
## Split an aligner SingleCellExperiment into one object per sampletag, using the
## demux_sampletags.R assignment. For each tag in use, writes the singlet cells
## (high-quality and called).
##
## --backend memory reads the source assays into memory once (sparse), then subsets in
## memory. --backend delayed keeps the assays on disk and subsets them lazily, which
## re-reads the source once per tag: low memory, but slow with many tags or large matrices.
##
## --output_format sce writes each split as an HDF5-backed SCE via
## saveHDF5SummarizedExperiment (se.rds marker; no top-level *_sce.rds, so the
## descriptive report's sce.rds glob skips them). h5ad writes each split as an h5ad via
## anndataR. both writes both. h5ad-only output cannot be checked by
## verify_sampletag_splits.R, which loads the HDF5 SCE.
##
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

script_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    f_flag <- which(ca == "-f")
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else if (length(f_flag) > 0 && length(ca) > f_flag[1]) {
        dirname(normalizePath(ca[f_flag[1] + 1L]))
    } else {
        "."
    }
})
source(file.path(script_dir, "sce_io.R"))

parser <- ArgumentParser(description = "Split an aligner SCE into per-sampletag SCEs.")
parser$add_argument("--sce", required = TRUE, help = "aligner SCE rds (HDF5-backed)")
parser$add_argument("--demux", required = TRUE, help = "demux assignment tsv.gz")
parser$add_argument("--out_base", required = TRUE, help = "base directory for by_sampletag output")
parser$add_argument("--backend", default = "memory", choices = c("memory", "delayed"),
                    help = "memory: realize assays in RAM once (fast); delayed: stream HDF5 per tag (low memory)")
parser$add_argument("--output_format", default = "sce", choices = c("sce", "h5ad", "both"),
                    help = "split file format(s) to write")
parser$add_argument("--tags", default = "", help = "optional comma-separated tag names to produce")
parser$add_argument("--labels", default = "", help = "optional comma-separated output labels, parallel to --tags")
parser$add_argument("--summary_rds", default = "", help = "optional summary rds for debugging")
parser$add_argument("--done", default = "", help = "optional sentinel file to touch when finished")
args <- parser$parse_args()

write_sce <- args$output_format %in% c("sce", "both")
write_h5ad <- args$output_format %in% c("h5ad", "both")

sce <- readRDS(args$sce)
if (args$backend == "memory") {
    sce <- realize_assays_in_memory(sce)
}

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
splits <- split_sce_singlets(sce, singlets, tags, labels)
summary <- data.frame(tag = character(), label = character(), n_cells = integer(),
                      stringsAsFactors = FALSE)

for (i in seq_along(tags)) {
    sub <- splits[[i]]
    out_dir <- file.path(args$out_base, labels[i])
    if (write_sce) saveHDF5SummarizedExperiment(sub, dir = out_dir, replace = TRUE)
    if (write_h5ad) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        write_sce_h5ad(sub, file.path(out_dir, "adata.h5ad"))
    }
    cat(sprintf("tag %s -> %s: %d cells\n", tags[i], labels[i], ncol(sub)))
    summary <- rbind(summary, data.frame(tag = tags[i], label = labels[i],
                                          n_cells = ncol(sub), stringsAsFactors = FALSE))
}

if (nzchar(args$summary_rds)) saveRDS(summary, args$summary_rds)
if (nzchar(args$done)) file.create(args$done)
