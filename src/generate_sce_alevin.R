#!/usr/bin/env R

suppressPackageStartupMessages( {
  library(SingleCellExperiment)
  library(argparse)
  library(Matrix)
  library(tximeta)
  library(DropletUtils)
})


parser <- ArgumentParser(description='Builds a WTA SingleCellExperiment object for a given sample - from alevin.')

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
                    help = 'native: tximeta filterBarcodes=TRUE; emptydrops: load unfiltered and apply DropletUtils emptyDrops')

args <- parser$parse_args()

wd <- args$working_dir
id <- args$sample

# read count matrix
dir<-file.path(wd, 'alevin', id)
files<-file.path(dir, "alevin", "quants_mat.gz")
file.exists(files)

filter_barcodes <- args$cell_filtering != 'emptydrops'
se <- tximeta(files, type="alevin", alevinArgs=list(filterBarcodes=filter_barcodes),
              txOut=TRUE, skipMeta=TRUE)

sce <- as(se, "SingleCellExperiment")

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
