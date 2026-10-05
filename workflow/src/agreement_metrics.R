## Agreement between two pseudobulk count vectors over the same genes.

## Mean absolute relative difference, mean(|x - y| / ((x + y) / 2)), over
## genes with a count in x or y. Symmetric in x and y.
mard <- function(x, y) {
    denom <- (x + y) / 2
    keep <- denom > 0
    mean(abs(x - y)[keep] / denom[keep], na.rm = TRUE)
}

counts_per_million <- function(x) x / sum(x) * 1e6

## MARD on counts per million, so a difference in total UMI yield cancels.
scaled_mard <- function(x, y) {
    mard(counts_per_million(x), counts_per_million(y))
}

spearman <- function(x, y) cor(x, y, method = "spearman")

## Spearman and scaled MARD of two named pseudobulk vectors on the genes both hold.
pseudobulk_agreement <- function(reference, other) {
    genes <- intersect(names(reference), names(other))
    data.frame(spearman = spearman(reference[genes], other[genes]),
               scaled_mard = scaled_mard(reference[genes], other[genes]))
}

## Square matrix of stat(mat[, a], mat[, b]) for every pair of columns.
pairwise_matrix <- function(mat, stat) {
    cols <- colnames(mat)
    out <- outer(cols, cols, Vectorize(function(a, b) stat(mat[, a], mat[, b])))
    dimnames(out) <- list(cols, cols)
    out
}

## One row per column pair, for plotting.
pairwise_long <- function(mat, stat) {
    long <- as.data.frame(as.table(pairwise_matrix(mat, stat)))
    names(long) <- c("pipeline1", "pipeline2", "value")
    long
}
