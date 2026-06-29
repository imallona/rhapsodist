#!/usr/bin/env R
##
## Check the per-sampletag splits from split_sce_by_sampletag.R. Loads each
## HDF5-backed SCE under a by_sampletag directory and stops on error unless every
## split is a valid SCE (HDF5-backed counts, same genes as the source) and the
## splits are disjoint with a union equal to the demux singlets present in the source.
##
## GPLv3

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(HDF5Array)
    library(data.table)
    library(argparse)
})

parser <- ArgumentParser(description = "Verify per-sampletag split SCE outputs.")
parser$add_argument("--by_dir", required = TRUE, help = "by_sampletag directory to check")
parser$add_argument("--sce", required = TRUE, help = "source aligner SCE rds")
parser$add_argument("--demux", required = TRUE, help = "demux assignment tsv.gz")
parser$add_argument("--expected_labels", default = "",
                    help = "optional comma-separated split dir names that must be present")
parser$add_argument("--out", default = "", help = "optional flag file to write on success")
args <- parser$parse_args()

source_sce <- readRDS(args$sce)
source_cells <- colnames(source_sce)
source_genes <- nrow(source_sce)

demux <- as.data.frame(fread(cmd = paste("zcat", shQuote(args$demux)), sep = "\t", header = TRUE))
singlets <- demux$cb[demux$status %in% c("highqual", "called") & !is.na(demux$called_tag)]
expected_cells <- intersect(singlets, source_cells)

split_dirs <- list.dirs(args$by_dir, recursive = FALSE, full.names = TRUE)
split_dirs <- split_dirs[file.exists(file.path(split_dirs, "se.rds"))]
if (length(split_dirs) == 0) {
    stop("no split objects found under ", args$by_dir)
}

if (nzchar(args$expected_labels)) {
    want <- strsplit(args$expected_labels, ",")[[1]]
    got <- basename(split_dirs)
    if (!setequal(want, got)) {
        stop("split dirs ", paste(sort(got), collapse = ", "),
             " do not match expected labels ", paste(sort(want), collapse = ", "))
    }
}

seen <- character(0)
total <- 0L
for (d in split_dirs) {
    sub <- loadHDF5SummarizedExperiment(d)
    label <- basename(d)
    if (!is(sub, "SingleCellExperiment")) stop(label, ": not a SingleCellExperiment")
    if (!"counts" %in% assayNames(sub)) stop(label, ": no counts assay")
    if (!is(assay(sub, "counts"), "DelayedArray")) stop(label, ": counts assay is not HDF5-backed")
    if (nrow(sub) != source_genes) {
        stop(label, ": ", nrow(sub), " genes, expected ", source_genes, " from the source SCE")
    }
    cells <- colnames(sub)
    if (!all(cells %in% source_cells)) stop(label, ": contains cells absent from the source SCE")
    if (length(intersect(cells, seen)) > 0) stop(label, ": shares cells with another split")
    seen <- c(seen, cells)
    total <- total + ncol(sub)
    cat(sprintf("%s: %d cells, %d genes, HDF5-backed counts\n", label, ncol(sub), nrow(sub)))
}

if (!setequal(seen, expected_cells)) {
    stop("union of split cells (", length(seen), ") does not match the demux singlet ",
         "cells present in the source SCE (", length(expected_cells), ")")
}
if (total == 0) stop("splits contain no cells")

cat(sprintf("OK: %d splits, %d cells total, partition matches the demux singlets\n",
            length(split_dirs), total))
if (nzchar(args$out)) file.create(args$out)
