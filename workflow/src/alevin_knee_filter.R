#!/usr/bin/env Rscript
## Filter alevin barcodes by a DropletUtils barcodeRanks knee on DeduplicatedReads.
## Usage: Rscript alevin_knee_filter.R \
##            --feature_dump <path> \
##            --output <path> \
##            [--sim_barcodes <path>]

suppressPackageStartupMessages({
    library(DropletUtils)
    library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

parse_args = function(argv) {
    out = list(feature_dump = NULL, output = NULL, sim_barcodes = "")
    i = 1
    while (i <= length(argv)) {
        if (argv[i] == "--feature_dump") { out$feature_dump = argv[i + 1]; i = i + 2 }
        else if (argv[i] == "--output")  { out$output = argv[i + 1];        i = i + 2 }
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

fd = fread(opts$feature_dump, sep = "\t", header = TRUE)
barcodes = fd[[1]]
dedup = as.numeric(fd[["DeduplicatedReads"]])

br = barcodeRanks(matrix(dedup, nrow = 1, dimnames = list(NULL, barcodes)))
knee_threshold = metadata(br)$knee

keep = dedup >= knee_threshold
writeLines(barcodes[keep], opts$output)

message(sprintf(
    "barcodeRanks knee: kept %d / %d barcodes (knee UMI threshold: %g)",
    sum(keep), length(barcodes), knee_threshold
))
