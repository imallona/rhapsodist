#!/usr/bin/env Rscript
## Self-contained tests for workflow/src/agreement_metrics.R.
## Run with:
##   cd workflow && Rscript tests/test_agreement_metrics.R

script_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
        "."
    }
})
source(file.path(script_dir, "..", "src", "agreement_metrics.R"))

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

near <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))

x <- c(10, 0, 4, 0, 7)
y <- c(5, 2, 4, 0, 21)

## per gene: 5/7.5, 2/1, 0/4, skipped, 14/14
check("mard matches the hand computation",
      near(mard(x, y), mean(c(5 / 7.5, 2, 0, 1))))
check("mard is symmetric", near(mard(x, y), mard(y, x)))
check("mard of a vector with itself is 0", near(mard(x, x), 0))
check("mard skips genes absent from both vectors",
      near(mard(c(1, 0), c(3, 0)), 1))
check("mard has an upper bound of 2", near(mard(c(1, 0), c(0, 1)), 2))

check("counts per million sum to 1e6", near(sum(counts_per_million(x)), 1e6))

check("a pure yield difference gives a non-zero mard", mard(x, 3 * x) > 0)
check("scaled mard is 0 for a pure yield difference",
      near(scaled_mard(x, 3 * x), 0))
check("scaled mard is symmetric", near(scaled_mard(x, y), scaled_mard(y, x)))

check("spearman is 1 for a monotone change", near(spearman(x, x^2 + 1), 1))
check("spearman is symmetric", near(spearman(x, y), spearman(y, x)))

mat <- cbind(a = x, b = y, c = 3 * x)
pm <- pairwise_matrix(mat, mard)
check("pairwise matrix keeps column names",
      identical(dimnames(pm), list(colnames(mat), colnames(mat))))
check("pairwise matrix is symmetric for mard", near(pm, t(pm)))
check("pairwise matrix has a zero diagonal for mard", near(unname(diag(pm)), c(0, 0, 0)))
check("pairwise matrix entry equals the direct call", near(pm["a", "b"], mard(x, y)))

long <- pairwise_long(mat, scaled_mard)
check("pairwise long has one row per column pair", nrow(long) == 9L)
check("pairwise long has the plotting columns",
      identical(names(long), c("pipeline1", "pipeline2", "value")))
check("pairwise long carries the statistic",
      near(long$value[long$pipeline1 == "a" & long$pipeline2 == "c"], 0))

cat("\nAll agreement_metrics tests passed.\n")
