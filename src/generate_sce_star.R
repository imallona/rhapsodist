#!/usr/bin/env Rscript
##
## Generates a SCE object with the relevant WTA data
##
## Izaskun Mallona
## Started Aug 12th 2024, reusing code from Oct 18th 2023

suppressPackageStartupMessages( {
    library(SingleCellExperiment)
    library(argparse)
    library(Matrix)
    library(DropletUtils)
})

parser <- ArgumentParser(description='Builds a WTA SingleCellExperiment object for a given sample.')

parser$add_argument('--sample',
                    type = "character",
                    help = 'Sample identifier')

## parser$add_argument('--run_mode', 
##                     type = 'character',
##                     help = 'Run mode')

parser$add_argument('--working_dir', 
                    type = 'character',
                    help = 'Working directory')

parser$add_argument('--output_fn',
                    type = 'character',
                    help = 'Output SCE filename (path)')

parser$add_argument('--cell_filtering',
                    type = 'character', default = 'native',
                    help = 'native: use pre-filtered STARsolo output; emptydrops: apply DropletUtils emptyDrops on raw counts')

## parser$add_argument('--captured_gtf', 
##                     type = 'character',
##                     help = 'Captured features GTF (path)')

args <- parser$parse_args()

## get_captured_gene_ids <- function(gtf) {
##     fd <- read.table(gtf, sep = '\t')
##     return(sapply(strsplit(fd$V9, split = ';'), function(x) return(gsub('gene_id ', '', x[1]))))
## }

read_matrix <- function(mtx, cells, features, cell.column = 1, feature.column = 1,
                        modality = 'wta') {
  cell.barcodes <- read.table(
    file = cells,
    header = FALSE,
    row.names = cell.column)


  feature.names <- read.table(
    file = features,
    header = FALSE,
    row.names = feature.column)

  d <- as(readMM(mtx), 'CsparseMatrix')

  ## if (modality == 'wta') {
  ## colnames(d) <- gsub('_', '', rownames(cell.barcodes))
  ## } else if (modality == 'tso') {
  ##     ## remove the fixed parts of the TSO CBs
  ##     colnames(d) <- paste0(
  ##         substr(rownames(cell.barcodes), 1, 9),
  ##         substr(rownames(cell.barcodes), 9+4+1, 9+4+9),
  ##         substr(rownames(cell.barcodes), 9+4+9+4+1, 9+4+9+4+9))
  ## }

  colnames(d) <- rownames(cell.barcodes)
  rownames(d) <- rownames(feature.names)
  
  return(d)
}

## wta start

wd <- args$working_dir
id <- args$sample

gene_dir <- if (args$cell_filtering == 'emptydrops') 'raw' else 'filtered'

wta <- read_matrix(mtx = file.path(wd, 'starsolo', id, 'Solo.out', 'Gene', gene_dir, 'matrix.mtx'),
                   cells = file.path(wd, 'starsolo', id, 'Solo.out', 'Gene', gene_dir, 'barcodes.tsv'),
                   features = file.path(wd, 'starsolo', id, 'Solo.out', 'Gene', gene_dir, 'features.tsv'),
                   cell.column = 1,
                   feature.column = 1)

wta_feat <- read.table(file.path(wd, 'starsolo', id, 'Solo.out', 'Gene', gene_dir, 'features.tsv'),
                       row.names = 1,
                       header = FALSE)

## captured <- get_captured_gene_ids(gtf)

colnames(wta_feat) <- c("name", "type", "value")

## wta_feat$captured <- ifelse(wta_feat$name %in% captured, yes = 'captured', no = 'not_captured')

## wta end

sce <- SingleCellExperiment(assays = list(counts = wta),
                            mainExpName = id,
                            rowData = wta_feat)

if (args$cell_filtering == 'emptydrops') {
    set.seed(42)
    ed <- tryCatch(
        emptyDrops(counts(sce)),
        error = function(e) {
            cat(sprintf('emptyDrops failed (%s); keeping all %d barcodes\n', conditionMessage(e), ncol(sce)))
            NULL
        }
    )
    if (!is.null(ed)) {
        keep <- !is.na(ed$FDR) & ed$FDR <= 0.01
        cat(sprintf('emptyDrops: kept %d / %d barcodes at FDR 0.01\n', sum(keep), ncol(sce)))
        sce <- sce[, keep]
    }
}

saveRDS(object = sce, file = args$output_fn)
