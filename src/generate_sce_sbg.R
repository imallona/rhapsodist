#!/usr/bin/env Rscript
##
## Generate a SingleCellExperiment object from SevenBridges / BD Rhapsody
## official pipeline MEX output.
##
## The official pipeline outputs cell barcodes as numeric indices (one-based,
## encoding the three-part CB1·CB2·CB3 combination).  This script decodes them
## to 27-bp concatenated sequences using the BD lookup tables in
## scripts/index2barcode.R so that barcodes are directly comparable to those
## produced by STARsolo, kallisto/bustools, and Alevin.
##
## Izaskun Mallona
## GPLv3

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(Matrix)
    library(argparse)
})

parser <- ArgumentParser(
    description = 'Build SCE from SBG/BD Rhapsody official pipeline MEX output.')

parser$add_argument('--sample',
    type = 'character', help = 'Sample identifier')
parser$add_argument('--mex_dir',
    type = 'character',
    help = 'Directory containing features.tsv[.gz], matrix.mtx[.gz], barcodes.tsv[.gz]')
parser$add_argument('--output_fn',
    type = 'character', help = 'Output RDS path')
parser$add_argument('--bead_version',
    type = 'character', default = 'Enh',
    help = 'BD bead version for barcode decoding: v1, Enh, or EnhV2 (default: Enh)')
parser$add_argument('--index2barcode_script',
    type = 'character',
    help = 'Path to src/index2barcode.R containing index_to_sequence()')
parser$add_argument('--whitelist_dir',
    type = 'character', default = NULL,
    help = 'Directory with BD_CLS1.txt, BD_CLS2.txt, BD_CLS3.txt (required for EnhV2)')
parser$add_argument('--features_map',
    type = 'character', default = NULL,
    help = 'STARsolo features.tsv (gene_id<TAB>gene_name<TAB>...) used to remap SBG gene symbols to Ensembl IDs')

args <- parser$parse_args()

if (identical(args$bead_version, 'EnhV2') && is.null(args$whitelist_dir)) {
    stop(
        "--whitelist_dir is required when --bead_version EnhV2. ",
        "Supply the directory containing BD_CLS1.txt, BD_CLS2.txt, BD_CLS3.txt."
    )
}

source(args$index2barcode_script)

## Load 384-sequence whitelist into B384_cell_key1/2/3 for EnhV2 decoding.
if (!is.null(args$whitelist_dir)) {
    B384_cell_key1 <<- readLines(file.path(args$whitelist_dir, 'BD_CLS1.txt'))
    B384_cell_key2 <<- readLines(file.path(args$whitelist_dir, 'BD_CLS2.txt'))
    B384_cell_key3 <<- readLines(file.path(args$whitelist_dir, 'BD_CLS3.txt'))
}

## locate MEX files, accepting both compressed and plain variants
.find <- function(dir, base) {
    gz  <- file.path(dir, paste0(base, '.gz'))
    raw <- file.path(dir, base)
    if (file.exists(gz))  return(gz)
    if (file.exists(raw)) return(raw)
    stop('Cannot find ', base, '[.gz] in ', dir)
}

mex_dir    <- args$mex_dir
feature_fn <- .find(mex_dir, 'features.tsv')
count_fn   <- .find(mex_dir, 'matrix.mtx')
barcode_fn <- .find(mex_dir, 'barcodes.tsv')

counts   <- readMM(count_fn)
features <- read.table(feature_fn, header = FALSE, sep = '\t',
                       stringsAsFactors = FALSE)
barcodes <- read.table(barcode_fn, header = FALSE, sep = '\t',
                       stringsAsFactors = FALSE)

rownames(counts) <- features$V1
colnames(counts) <- barcodes$V1

sce <- SingleCellExperiment(assays = list(counts = counts),
                            mainExpName = args$sample)

## Decode numeric barcode indices → 27-bp concatenated CB1·CB2·CB3 sequences
colnames(sce) <- sapply(colnames(sce), index_to_sequence,
                        bead_version = args$bead_version)

## Remap SBG gene symbols to Ensembl IDs using STARsolo features.tsv so that
## rownames match the other aligners in the comparison report.
if (!is.null(args$features_map) && file.exists(args$features_map)) {
    feat <- read.table(args$features_map, header = FALSE, sep = '\t',
                       stringsAsFactors = FALSE)
    ## col 1 = gene_id (Ensembl), col 2 = gene_name (symbol)
    name_to_id <- setNames(feat$V1, feat$V2)
    remapped <- name_to_id[rownames(sce)]
    valid <- !is.na(remapped)
    cat(sprintf('Gene symbol remapping: %d / %d genes matched to Ensembl IDs\n',
                sum(valid), nrow(sce)))
    rownames(sce)[valid] <- remapped[valid]
} else {
    cat('No features_map provided; keeping gene symbols as rownames\n')
}

dir.create(dirname(args$output_fn), recursive = TRUE, showWarnings = FALSE)
saveRDS(sce, args$output_fn)
