#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(argparse)
    library(Matrix)
    library(DropletUtils)
    library(HDF5Array)
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
source(file.path(script_dir, "gtf_utils.R"))

parser <- ArgumentParser(description='Builds a WTA SingleCellExperiment object for a given sample - from Kallisto.')

parser$add_argument('--sample',
                    type = "character",
                    help = 'Sample identifier')

parser$add_argument('--working_dir',
                    type = 'character',
                    help = 'Working directory')

parser$add_argument('--gtf',
                    type = 'character',
                    help = 'GTF annotation used to map gene_id to gene_name')

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

gtf_genes <- parse_gtf_genes(args$gtf)
stopifnot(nrow(gtf_genes) > 0)
row_data <- build_rowdata_from_gtf(gene_ids, gtf_genes)
stopifnot(identical(rownames(row_data), gene_ids),
          all(c("name", "type", "value") %in% colnames(row_data)))
matched <- sum(row_data$name != gene_ids)
cat(sprintf('kallisto gene symbol mapping: %d / %d gene_ids matched a GTF gene_name\n',
            matched, length(gene_ids)))
if (matched == 0L) {
    warning('no kallisto gene_ids matched any GTF gene_name; check --gtf and ID version suffixes')
}

# bustools output: rows = barcodes, cols = genes; SCE convention: rows = genes
sce <- SingleCellExperiment(list(counts = t(counts)),
                            colData = DataFrame(Barcode = barcodes),
                            rowData = row_data,
                            mainExpName = id)
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
