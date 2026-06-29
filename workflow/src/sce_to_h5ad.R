#!/usr/bin/env R
##
## Convert a main aligner SingleCellExperiment (rds, HDF5-backed) to an h5ad file via
## anndataR, for python-side or anndataR consumers. Additive: the source rds and its
## _hdf5 directory are left untouched.
##
## GPLv3

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(HDF5Array)
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

parser <- ArgumentParser(description = "Convert an aligner SCE rds to h5ad via anndataR.")
parser$add_argument("--sce", required = TRUE, help = "aligner SCE rds (HDF5-backed)")
parser$add_argument("--out", required = TRUE, help = "output h5ad path")
args <- parser$parse_args()

sce <- readRDS(args$sce)
write_sce_h5ad(sce, args$out)
cat(sprintf("wrote %s (%d genes x %d cells)\n", args$out, nrow(sce), ncol(sce)))
