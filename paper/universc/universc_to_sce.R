#!/usr/bin/env Rscript
## Writes the filtered UniverSC (Cell Ranger) matrix as an HDF5-backed
## SingleCellExperiment, laid out like the rhapsodist aligner outputs: rownames
## are gene IDs, rowData$name the gene symbol, colnames the 27 nt barcodes
## without Cell Ranger's "-1" suffix.
##
## Rscript paper/universc/universc_to_sce.R OUTS_DIR SAMPLE OUTPUT_RDS

suppressPackageStartupMessages({
    library(DropletUtils)
    library(HDF5Array)
    library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
outs_dir <- args[1]
sample <- args[2]
output_rds <- args[3]

sce <- read10xCounts(file.path(outs_dir, "filtered_feature_bc_matrix"),
                     col.names = TRUE)
rownames(sce) <- rowData(sce)$ID
rowData(sce) <- S4Vectors::DataFrame(name = rowData(sce)$Symbol,
                                     type = rowData(sce)$Type,
                                     row.names = rownames(sce))
colData(sce) <- NULL
colnames(sce) <- sub("-1$", "", colnames(sce))
mainExpName(sce) <- sample

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
sce <- saveHDF5SummarizedExperiment(sce, dir = sub("[.]rds$", "_hdf5", output_rds),
                                    replace = TRUE)
base::saveRDS(sce, output_rds)
