#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(argparse)
    library(Matrix)
    library(DropletUtils)
    library(HDF5Array)
})

parser <- ArgumentParser(description='Builds a WTA SingleCellExperiment object for a given sample - from Kallisto.')

parser$add_argument('--sample',
                    type = "character",
                    help = 'Sample identifier')

parser$add_argument('--working_dir', 
                    type = 'character',
                    help = 'Working directory')

parser$add_argument('--output_fn',
                    type = 'character',
                    help = 'Output SCE filename (path)')

parser$add_argument('--cell_filtering',
                    type = 'character', default = 'native',
                    help = 'native: barcodeRanks knee filter; emptydrops: apply DropletUtils emptyDrops')

args <- parser$parse_args()


wd <- args$working_dir
id <- args$sample

counts <- Matrix::readMM(file.path(wd, 'bustools', id, 'output.mtx'))
gene_ids <- readLines(file.path(wd, 'bustools', id, 'output.genes.txt'))
barcodes <- readLines(file.path(wd, 'bustools', id, 'output.barcodes.txt'))

# bustools output: rows = barcodes, cols = genes; SCE convention: rows = genes
sce <- SingleCellExperiment(list(counts = t(counts)),
                            colData = DataFrame(Barcode = barcodes),
                            rowData = DataFrame(ID = gene_ids, SYMBOL = gene_ids))
rownames(sce) <- gene_ids
colnames(sce) <- barcodes

saveRDS(sce, file.path(dirname(args$output_fn), paste0(id, '_kallisto_sce_pre_filter.rds')))

if (args$cell_filtering == 'native') {
    br <- barcodeRanks(counts(sce))
    knee_threshold <- metadata(br)$knee
    keep <- colSums(counts(sce)) >= knee_threshold
    cat(sprintf('barcodeRanks knee: kept %d / %d barcodes (knee UMI threshold: %g)\n',
                sum(keep), ncol(sce), knee_threshold))
    sce <- sce[, keep]
} else if (args$cell_filtering == 'emptydrops') {
    set.seed(42)
    ed <- tryCatch(
        emptyDrops(counts(sce)),
        error = function(e) {
            warning(sprintf('emptyDrops failed (%s); keeping all %d barcodes', conditionMessage(e), ncol(sce)))
            NULL
        }
    )
    if (!is.null(ed)) {
        keep <- !is.na(ed$FDR) & ed$FDR <= 0.01
        cat(sprintf('emptyDrops: kept %d / %d barcodes at FDR 0.01\n', sum(keep), ncol(sce)))
        sce <- sce[, keep]
    }
}

hdf5_dir <- sub('\\.rds$', '_hdf5', args$output_fn)
sce <- saveHDF5SummarizedExperiment(sce, dir = hdf5_dir, replace = TRUE)
base::saveRDS(sce, args$output_fn)
