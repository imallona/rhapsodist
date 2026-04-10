#!/usr/bin/env Rscript
## Filter alevin barcodes using DropletUtils barcodeRanks.
## Non-sketch mode: reads DeduplicatedReads from featureDump.txt.
## Sketch mode: reads per-cell UMI totals from the alevin-fry MTX output.
## Usage: Rscript alevin_knee_filter.R \
##            (--feature_dump <path> | --fry_quant_dir <path>) \
##            --output <path> \
##            [--sim_barcodes <path>]

suppressPackageStartupMessages({
    library(DropletUtils)
    library(Matrix)
    library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

parse_args = function(argv) {
    out = list(feature_dump = "", fry_quant_dir = "", output = NULL, sim_barcodes = "")
    i = 1
    while (i <= length(argv)) {
        if (argv[i] == "--feature_dump")   { out$feature_dump = argv[i + 1];   i = i + 2 }
        else if (argv[i] == "--fry_quant_dir") { out$fry_quant_dir = argv[i + 1]; i = i + 2 }
        else if (argv[i] == "--output")    { out$output = argv[i + 1];         i = i + 2 }
        else if (argv[i] == "--sim_barcodes") { out$sim_barcodes = argv[i + 1]; i = i + 2 }
        else i = i + 1
    }
    out
}

opts = parse_args(args)

if (nchar(opts$sim_barcodes) > 0) {
    file.copy(opts$sim_barcodes, opts$output, overwrite = TRUE)
    message("simulation mode: copied ", opts$sim_barcodes, " to ", opts$output)
    quit(save = "no", status = 0)
}

if (nchar(opts$fry_quant_dir) > 0) {
    fry_alevin = file.path(opts$fry_quant_dir, "alevin")
    mtx_f = file.path(fry_alevin, "quants_mat.mtx")
    rows_f = file.path(fry_alevin, "quants_mat_rows.txt")

    mat = readMM(mtx_f)
    barcodes = readLines(rows_f)
    umi_totals = rowSums(mat)
} else {
    fd = fread(opts$feature_dump, sep = "\t", header = TRUE)
    barcodes = fd[[1]]
    umi_totals = as.numeric(fd[["DeduplicatedReads"]])
}

br = barcodeRanks(matrix(umi_totals, nrow = 1, dimnames = list(NULL, barcodes)))
knee_threshold = metadata(br)$knee

keep = umi_totals >= knee_threshold
writeLines(barcodes[keep], opts$output)

message(sprintf(
    "barcodeRanks knee: kept %d / %d barcodes (knee UMI threshold: %g)",
    sum(keep), length(barcodes), knee_threshold
))
