#!/usr/bin/env Rscript
## Tests for paper/compare_workflows.R. Run from the repository root:
##   Rscript paper/test_compare_workflows.R

sys.source("paper/compare_workflows.R", envir = environment())

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

sparse <- function(values, genes, cells) {
    Matrix::Matrix(values, nrow = length(genes), sparse = TRUE, dimnames = list(genes, cells))
}
aligner_summaries <- list(
    starsolo = tool_summary(sparse(c(1, 2, 3, 4), c("g1", "g2"), c("a", "b"))),
    alevin = tool_summary(sparse(c(4, 3, 2, 1, 0, 0), c("g1", "g2"), c("a", "b", "c"))))
benchmark <- data.table(seconds = 60, max_rss_mb = 100)

other <- tool_summary(sparse(c(2, 4, 6, 8, 1, 1), c("g1", "g2"), c("a", "b", "z")))
check("a summary holds barcodes and gene totals",
      identical(other$cells, c("a", "b", "z")) && identical(other$pseudobulk, c(g1 = 9, g2 = 13)))
rows <- workflow_rows("other", other, aligner_summaries, benchmark)
check("one row per aligner", identical(rows$aligner, c("starsolo", "alevin")))
check("cells and status of a finished workflow",
      all(rows$status == "ok") && all(rows$n_cells == 3))
check("cells shared by barcode", identical(rows$shared_cells, c(2L, 2L)))
check("agreement is computed per aligner",
      isTRUE(all.equal(rows$spearman, c(1, -1))))
check("time and memory are repeated on every row",
      all(rows$seconds == 60) && all(rows$max_rss_mb == 100))

failed <- workflow_rows("other", NULL, aligner_summaries, benchmark)
check("a workflow without counts is listed as failed",
      all(failed$status == "failed") && all(is.na(failed$n_cells)) && all(is.na(failed$spearman)))

absent <- read_benchmark(file.path(tempdir(), "no_such_benchmark.txt"))
check("a missing benchmark file gives NA", is.na(absent$seconds) && is.na(absent$max_rss_mb))
path <- tempfile()
fwrite(data.table(s = c(10, 5), max_rss = c(7, 9)), path, sep = "\t")
present <- read_benchmark(path)
check("seconds are summed and memory is the peak",
      present$seconds == 15 && present$max_rss_mb == 9)
