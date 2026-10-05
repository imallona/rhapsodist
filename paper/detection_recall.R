## Detection of gene-by-cell entries against a known count matrix.

suppressPackageStartupMessages(library(data.table))

TRUE_COUNT_BREAKS <- c(-Inf, 0, 1, 2, 5, 10, Inf)
TRUE_COUNT_LABELS <- c("0", "1", "2", "3-5", "6-10", ">10")

## Fraction of entries with an observed count above zero, per bin of the true
## count, on the genes and cells both matrices hold. For the bin "0" this is
## false detection; for the other bins it is recall.
detection_by_true_count <- function(truth, observed) {
    genes <- intersect(rownames(truth), rownames(observed))
    cells <- intersect(colnames(truth), colnames(observed))
    true_counts <- as.vector(as.matrix(truth[genes, cells, drop = FALSE]))
    detected <- as.vector(as.matrix(observed[genes, cells, drop = FALSE])) > 0
    bins <- cut(true_counts, breaks = TRUE_COUNT_BREAKS, labels = TRUE_COUNT_LABELS)
    data.table(true_count = factor(TRUE_COUNT_LABELS, levels = TRUE_COUNT_LABELS),
               n_entries = as.vector(table(bins)),
               detected_fraction = as.vector(tapply(detected, bins, mean)))
}
