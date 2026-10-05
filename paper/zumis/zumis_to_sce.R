#!/usr/bin/env Rscript
## Writes the zUMIs exon UMI counts of the barcodes zUMIs kept as an
## HDF5-backed SingleCellExperiment, laid out like the rhapsodist aligner
## outputs: rownames are gene IDs, rowData$name the gene symbol.
##
## Rscript paper/zumis/zumis_to_sce.R DGECOUNTS_RDS SAMPLE OUTPUT_RDS

suppressPackageStartupMessages({
    library(HDF5Array)
    library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
counts_rds <- args[1]
sample <- args[2]
output_rds <- args[3]

counts <- readRDS(counts_rds)$umicount$exon$all
gene_names <- read.delim(sub("dgecounts[.]rds$", "gene_names.txt", counts_rds))
symbols <- gene_names$gene_name[match(rownames(counts), gene_names$gene_id)]

sce <- SingleCellExperiment(assays = list(counts = counts),
                            rowData = S4Vectors::DataFrame(name = symbols,
                                                           row.names = rownames(counts)))
mainExpName(sce) <- sample

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
sce <- saveHDF5SummarizedExperiment(sce, dir = sub("[.]rds$", "_hdf5", output_rds),
                                    replace = TRUE)
base::saveRDS(sce, output_rds)
