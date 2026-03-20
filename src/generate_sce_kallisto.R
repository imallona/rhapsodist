#!/usr/bin/env Rscript

suppressPackageStartupMessages( {
  library(SingleCellExperiment)
  library(argparse)
  library(Matrix)
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

rownames(counts) <- barcodes
colnames(counts) <- gene_ids

# Cell filtering: elbow detection in log-log barcode rank space.
# A min-count pre-filter (>= 5 UMIs) is applied before elbow detection to
# ignore the massive tail of likely-empty barcodes from unfiltered bustools
# output, then the maximum-distance-from-diagonal method finds the knee.
lib_sizes <- Matrix::rowSums(counts)
nonzero <- lib_sizes[lib_sizes > 0]
if (length(nonzero) == 0L) {
    stop("No non-zero library sizes detected; kallisto output appears empty or invalid.")
}
candidates <- nonzero[nonzero >= 5]
if (length(candidates) < 10L) {
    candidates <- nonzero
}
ranked <- sort(candidates, decreasing = TRUE)
n <- length(ranked)
log_rank <- log10(seq_len(n))
log_count <- log10(ranked)
x1 <- log_rank[1]; y1 <- log_count[1]
x2 <- log_rank[n]; y2 <- log_count[n]
dx <- x2 - x1; dy <- y2 - y1
dist <- (dy * log_rank - dx * log_count + x2*y1 - y2*x1) / sqrt(dy^2 + dx^2)
knee_idx <- which.max(dist)
knee_threshold <- ranked[knee_idx]
n_before <- length(lib_sizes)
keep <- lib_sizes >= knee_threshold
cat(sprintf('Elbow threshold: %g counts  Kept %d / %d barcodes\n',
            knee_threshold, sum(keep), n_before))

counts <- counts[keep, ]
barcodes <- barcodes[keep]

sce <- SingleCellExperiment(list(counts = t(counts)),
                            colData = DataFrame(Barcode = barcodes),
                            rowData = DataFrame(ID = gene_ids, SYMBOL = gene_ids))
rownames(sce) <- gene_ids

saveRDS(object = sce, file = args$output_fn)
