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

# Cell filtering via DropletUtils::barcodeRanks().
# bustools count produces unfiltered output (no bustools correct step), so
# the matrix contains all barcodes including empty droplets.  barcodeRanks
# fits a smooth rank-count curve and returns the inflection point, which is
# more robust on the large unfiltered pools than a simple diagonal heuristic.
br <- barcodeRanks(counts(sce))
knee_threshold <- metadata(br)$inflection
n_before <- ncol(sce)
keep <- colSums(counts(sce)) >= knee_threshold
cat(sprintf('DropletUtils inflection threshold: %g counts  Kept %d / %d barcodes\n',
            knee_threshold, sum(keep), n_before))
sce <- sce[, keep]

saveRDS(object = sce, file = args$output_fn)
