#!/usr/bin/env R
##
## Shared SCE input/output helpers for the sampletag split and the h5ad export.
## No side effects at source time, so the functions can be unit tested.
##
## GPLv3

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
})

## Coerce every assay to an in-memory CsparseMatrix. For an HDF5-backed SCE this reads
## the source once, instead of leaving DelayedArrays that get re-read on each later
## subset. This is the 'memory' split backend.
realize_assays_in_memory <- function(sce) {
    for (a in assayNames(sce)) {
        assay(sce, a, withDimnames = FALSE) <- as(assay(sce, a, withDimnames = FALSE),
                                                  "CsparseMatrix")
    }
    sce
}

## Partition an SCE into one subset per tag using the demux singlet assignment.
## singlets is a data frame with columns cb (cell barcode) and called_tag. Returns a
## named list (names are labels) of column subsets restricted to cells present in sce.
## Pure logic: no realization or I/O, so the partition is testable without HDF5.
split_sce_singlets <- function(sce, singlets, tags, labels) {
    stopifnot(length(tags) == length(labels))
    present <- colnames(sce)
    splits <- vector("list", length(tags))
    names(splits) <- labels
    for (i in seq_along(tags)) {
        barcodes <- singlets$cb[singlets$called_tag == tags[i]]
        cells <- intersect(barcodes, present)
        splits[[i]] <- sce[, cells, drop = FALSE]
    }
    splits
}

## Write an SCE to an h5ad file via anndataR, which needs in-memory matrices, so the
## assays are realized first.
write_sce_h5ad <- function(sce, path) {
    if (!requireNamespace("anndataR", quietly = TRUE)) {
        stop("anndataR is required to write h5ad output; add it to the r_bioc environment")
    }
    adata <- anndataR::as_AnnData(realize_assays_in_memory(sce))
    anndataR::write_h5ad(adata, path, mode = "w")
    invisible(path)
}
