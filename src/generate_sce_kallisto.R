#!/usr/bin/env Rscript

suppressPackageStartupMessages( {
  library(SingleCellExperiment)
  library(argparse)
  library(Matrix)
  library(DropletUtils)
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
                    help = 'native: no filtering applied (raw bustools output); emptydrops: apply DropletUtils emptyDrops')

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

if (args$cell_filtering == 'emptydrops') {
    set.seed(42)
    ed <- emptyDrops(counts(sce))
    keep <- !is.na(ed$FDR) & ed$FDR <= 0.01
    cat(sprintf('emptyDrops: kept %d / %d barcodes at FDR 0.01\n', sum(keep), ncol(sce)))
    sce <- sce[, keep]
}

saveRDS(object = sce, file = args$output_fn)
