#!/usr/bin/env R

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(argparse)
    library(Matrix)
    library(DropletUtils)
    library(HDF5Array)
})

parser <- ArgumentParser(
    description = 'Build a WTA SingleCellExperiment from alevin output, using a pre-computed knee barcode list.')

parser$add_argument('--sample', type = 'character', help = 'Sample identifier')
parser$add_argument('--working_dir', type = 'character', help = 'Working directory')
parser$add_argument('--output_fn', type = 'character', help = 'Output SCE RDS path')
parser$add_argument('--knee_barcodes', type = 'character',
                    help = 'Path to knee-filtered barcode list from alevin_knee_filter.py')
parser$add_argument('--cell_filtering', type = 'character', default = 'native',
                    help = 'native: use knee barcodes directly; emptydrops: apply DropletUtils emptyDrops')

args <- parser$parse_args()

wd <- args$working_dir
id <- args$sample
alevin_dir <- file.path(wd, 'alevin', id, 'alevin')

rows_file <- file.path(alevin_dir, 'quants_mat_rows.txt')
cols_file <- file.path(alevin_dir, 'quants_mat_cols.txt')
matrix_file <- file.path(alevin_dir, 'quants_mat.gz')

all_barcodes <- readLines(rows_file)
gene_names <- readLines(cols_file)
knee_barcodes <- readLines(args$knee_barcodes)

# streaming binary reader: reads quants_mat.gz row by row, keeping only knee barcodes.
# the alevin binary format stores per barcode: a bit-packed presence vector (ceil(n_genes/8) bytes)
# followed by float32 counts for each expressed gene, in order.
read_alevin_filtered <- function(matrix_file, gene_names, all_barcodes, keep_barcodes) {
    n_genes <- length(gene_names)
    n_all <- length(all_barcodes)
    len_bit_vec <- ceiling(n_genes / 8L)
    keep_set <- keep_barcodes[keep_barcodes %in% all_barcodes]
    keep_lookup <- all_barcodes %in% keep_set
    n_keep <- sum(keep_lookup)

    gene_idx_list <- vector('list', n_keep)
    count_list <- vector('list', n_keep)
    kept_names <- character(n_keep)
    pos <- 0L

    con <- gzcon(file(matrix_file, 'rb'))
    for (j in seq_len(n_all)) {
        ints <- readBin(con, integer(), size = 1L, signed = FALSE,
                        endian = 'little', n = len_bit_vec)
        bits <- matrix(intToBits(ints), nrow = 32L)
        mode(bits) <- 'integer'
        bits <- bits[8:1, ]
        expressed <- which(head(as.vector(bits), n_genes) == 1L)
        n_exp <- length(expressed)
        cts <- readBin(con, double(), size = 4L, endian = 'little', n = n_exp)
        if (keep_lookup[j]) {
            pos <- pos + 1L
            gene_idx_list[[pos]] <- expressed
            count_list[[pos]] <- cts
            kept_names[pos] <- all_barcodes[j]
        }
    }
    close(con)

    lens <- lengths(gene_idx_list)
    sparseMatrix(
        i = unlist(gene_idx_list),
        j = rep(seq_len(n_keep), lens),
        x = unlist(count_list),
        dims = c(n_genes, n_keep),
        dimnames = list(gene_names, kept_names),
        repr = 'T'
    )
}

cat('reading alevin binary matrix for', length(knee_barcodes), 'knee barcodes out of',
    length(all_barcodes), 'total\n')
mat <- read_alevin_filtered(matrix_file, gene_names, all_barcodes, knee_barcodes)

sce <- SingleCellExperiment(list(counts = mat), mainExpName = id)

saveRDS(sce, file.path(dirname(args$output_fn), paste0(id, '_alevin_sce_pre_filter.rds')))

if (args$cell_filtering == 'emptydrops') {
    set.seed(42)
    ed <- tryCatch(
        emptyDrops(counts(sce)),
        error = function(e) {
            cat(sprintf('emptyDrops failed (%s); keeping all %d barcodes\n',
                        conditionMessage(e), ncol(sce)))
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
