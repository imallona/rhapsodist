#!/usr/bin/env Rscript
## Compares rhapsodist sample tag calls with per-tag expression files of the BD
## pipeline (columns Cell_Index, Bioproduct, ...), one file per tag.
##
## Rscript paper/compare_sampletag_calls.R \
##     --demux output/gse282765/sampletags/colon/sampletag_demux.tsv.gz \
##     --deposited mouse_sampletag_3=PHIL_ST03.st.txt.gz mouse_sampletag_4=WT_ST04.st.txt.gz \
##     --whitelist_dir workflow/data/whitelist_384x3 \
##     --out_prefix output/gse282765/paper/colon/sampletag_concordance

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(ggplot2)
})

## BD cell index to the 27 nt barcode of the three whitelist segments.
cell_index_to_barcode <- function(index, whitelists) {
    n <- length(whitelists[[1]])
    zero_based <- index - 1
    paste0(whitelists[[1]][(zero_based %/% n %/% n) %% n + 1],
           whitelists[[2]][(zero_based %/% n) %% n + 1],
           whitelists[[3]][zero_based %% n + 1])
}

read_whitelists <- function(whitelist_dir) {
    lapply(1:3, function(i) readLines(file.path(whitelist_dir, sprintf("BD_CLS%d.txt", i))))
}

read_deposited_cells <- function(path) {
    indices <- fread(cmd = paste("zcat", shQuote(path), "| grep -v '^#' | cut -f1 | uniq"),
                     header = TRUE)[[1]]
    unique(indices)
}

## One row per deposited cell: its deposited tag and the rhapsodist outcome
## (same tag, other tag, multiplet, undetermined, or absent from the demux table).
classify_calls <- function(deposited, demux) {
    calls <- merge(deposited, demux[, .(cb, called_tag, status)], by = "cb", all.x = TRUE)
    calls[, outcome := fifelse(is.na(status), "absent",
                       fifelse(status %in% c("multiplet", "undetermined"), status,
                       fifelse(called_tag == deposited_tag, "same tag", "other tag")))]
    calls[]
}

main <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--demux", required = TRUE, help = "sampletag_demux.tsv.gz of a sample")
    parser$add_argument("--deposited", required = TRUE, nargs = "+",
                        help = "tag=file pairs, one BD expression file per tag")
    parser$add_argument("--whitelist_dir", required = TRUE)
    parser$add_argument("--out_prefix", required = TRUE)
    args <- parser$parse_args()

    whitelists <- read_whitelists(args$whitelist_dir)
    pairs <- strsplit(args$deposited, "=", fixed = TRUE)
    deposited <- rbindlist(lapply(pairs, function(pair) {
        data.table(deposited_tag = pair[1],
                   cb = cell_index_to_barcode(read_deposited_cells(pair[2]), whitelists))
    }))
    calls <- classify_calls(deposited, fread(cmd = paste("zcat", shQuote(args$demux))))

    counts <- calls[, .N, by = .(deposited_tag, outcome)][order(deposited_tag, -N)]
    counts[, fraction := round(N / sum(N), 4), by = deposited_tag]
    overall <- calls[, .(deposited_tag = "all", N = .N), by = outcome]
    overall[, fraction := round(N / sum(N), 4)]
    counts <- rbind(counts, overall[order(-N)], use.names = TRUE)

    dir.create(dirname(args$out_prefix), recursive = TRUE, showWarnings = FALSE)
    fwrite(counts, paste0(args$out_prefix, ".csv"))
    plot <- ggplot(counts[deposited_tag != "all"], aes(outcome, deposited_tag, fill = fraction)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = N), size = 3) +
        scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1)) +
        theme_bw() +
        labs(x = "rhapsodist call", y = "deposited tag", fill = "fraction")
    ggsave(paste0(args$out_prefix, ".pdf"), plot, width = 5, height = 2.5)
    print(counts)
}

if (sys.nframe() == 0) main()
