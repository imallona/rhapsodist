#!/usr/bin/env R

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

parser <- ArgumentParser(
    description = 'Build a WTA SingleCellExperiment from alevin output.')

parser$add_argument('--sample', type = 'character', help = 'Sample identifier')
parser$add_argument('--working_dir', type = 'character', default = '',
                    help = 'Working directory (non-sketch mode)')
parser$add_argument('--fry_quant_dir', type = 'character', default = '',
                    help = 'alevin-fry quant output directory (sketch mode)')
parser$add_argument('--gtf', type = 'character',
                    help = 'GTF annotation used to map gene_id to gene_name')
parser$add_argument('--output_fn', type = 'character', help = 'Output SCE RDS path')
parser$add_argument('--knee_barcodes', type = 'character',
                    help = 'Path to knee-filtered barcode list')
parser$add_argument('--cell_filtering', type = 'character', default = 'native',
                    help = 'native: use knee barcodes directly; emptydrops: apply DropletUtils emptyDrops; none: keep all barcodes')
parser$add_argument('--usa', action = 'store_true', default = FALSE,
                    help = 'alevin-fry USA mode: split spliced/unspliced/ambiguous counts into separate assays')

args <- parser$parse_args()

id <- args$sample
knee_barcodes <- readLines(args$knee_barcodes)

if (nchar(args$fry_quant_dir) > 0) {
    fry_alevin = file.path(args$fry_quant_dir, "alevin")
    mat_raw = readMM(file.path(fry_alevin, "quants_mat.mtx"))
    all_barcodes = readLines(file.path(fry_alevin, "quants_mat_rows.txt"))
    gene_names = readLines(file.path(fry_alevin, "quants_mat_cols.txt"))
    rownames(mat_raw) = all_barcodes
    colnames(mat_raw) = gene_names
    keep = all_barcodes %in% knee_barcodes
    cat('alevin-fry MTX: keeping', sum(keep), 'of', length(all_barcodes), 'barcodes\n')
    mat_raw = mat_raw[keep, , drop = FALSE]
    mat = t(mat_raw)
    colnames(mat) = rownames(mat_raw)
    rownames(mat) = gene_names
} else {
    wd <- args$working_dir
    alevin_dir <- file.path(wd, 'alevin', id, 'alevin')

    all_barcodes <- readLines(file.path(alevin_dir, 'quants_mat_rows.txt'))
    gene_names <- readLines(file.path(alevin_dir, 'quants_mat_cols.txt'))
    matrix_file <- file.path(alevin_dir, 'quants_mat.gz')

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
}

## USA mode: alevin-fry names the gene columns by splicing status: the spliced
## column is the bare gene_id, the unspliced column is <gene_id>-U and the
## ambiguous column is <gene_id>-A (there is no -S suffix). Collapse to one row
## per gene with separate assays: the main 'counts' is spliced plus ambiguous
## (the standard gene count), and spliced, unspliced and ambiguous are kept for
## RNA velocity.
collapse_usa_assays <- function(mat) {
    rn <- rownames(mat)
    is_unspliced <- grepl('-U$', rn)
    is_ambiguous <- grepl('-A$', rn)
    is_spliced <- !is_unspliced & !is_ambiguous
    base_id <- sub('-[UA]$', '', rn)
    genes <- sort(unique(base_id))
    cells <- colnames(mat)
    per_status <- function(sel) {
        out <- Matrix(0, nrow = length(genes), ncol = length(cells), sparse = TRUE,
                      dimnames = list(genes, cells))
        if (any(sel)) {
            out[match(base_id[sel], genes), ] <- mat[sel, , drop = FALSE]
        }
        out
    }
    spliced <- per_status(is_spliced)
    unspliced <- per_status(is_unspliced)
    ambiguous <- per_status(is_ambiguous)
    list(counts = spliced + ambiguous, spliced = spliced,
         unspliced = unspliced, ambiguous = ambiguous)
}

if (isTRUE(args$usa)) {
    assays_list <- collapse_usa_assays(mat)
    mat <- assays_list$counts
    cat(sprintf('USA mode: collapsed to %d genes; assays counts(=S+A), spliced, unspliced, ambiguous\n',
                nrow(mat)))
} else {
    assays_list <- list(counts = mat)
}

alevin_gene_ids <- rownames(mat)
gtf_genes <- parse_gtf_genes(args$gtf)
stopifnot(nrow(gtf_genes) > 0)
row_data <- build_rowdata_from_gtf(alevin_gene_ids, gtf_genes)
stopifnot(identical(rownames(row_data), alevin_gene_ids),
          all(c("name", "type", "value") %in% colnames(row_data)))
matched <- sum(row_data$name != alevin_gene_ids)
cat(sprintf('alevin gene symbol mapping: %d / %d gene_ids matched a GTF gene_name\n',
            matched, length(alevin_gene_ids)))
if (matched == 0L) {
    warning('no alevin gene_ids matched any GTF gene_name; check --gtf and ID version suffixes')
}

sce <- SingleCellExperiment(assays_list,
                            rowData = row_data,
                            mainExpName = id)

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
