#!/usr/bin/env Rscript
## Tests for paper/compare_linker_tolerance.R. Run from the repository root:
##   Rscript paper/test_compare_linker_tolerance.R

sys.source("paper/compare_linker_tolerance.R", envir = environment())

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

counts <- Matrix::Matrix(c(1, 0, 3,
                           2, 2, 0), nrow = 2, byrow = TRUE, sparse = TRUE,
                         dimnames = list(c("g1", "g2"), c("a", "b", "c")))
summary_row <- run_summary(counts)
check("cells are the matrix columns", summary_row$n_cells == 3)
check("median UMIs per cell", summary_row$median_umis == 3)

reference <- c(g1 = 10, g2 = 20, g3 = 30)
check("a run agrees with itself",
      identical(unlist(pseudobulk_agreement(reference, reference)),
                c(spearman = 1, scaled_mard = 0)))
doubled <- pseudobulk_agreement(reference, c(g3 = 60, g2 = 40, g1 = 20, g4 = 5))
check("genes are matched by name and a change in yield alone leaves scaled MARD at 0",
      doubled$spearman == 1 && isTRUE(all.equal(doubled$scaled_mard, 0)))
changed <- pseudobulk_agreement(reference, c(g1 = 30, g2 = 20, g3 = 10))
check("a reversed order gives Spearman -1 and a positive scaled MARD",
      changed$spearman == -1 && changed$scaled_mard > 0)

row_data <- S4Vectors::DataFrame(name = c("A", "B"))
summed <- symbol_pseudobulk(counts, row_data)
check("pseudobulk is named by gene symbol and sums the called cells",
      identical(summed, c(A = 4, B = 4)))
empty_cell <- cbind(counts, d = c(0, 0))
check("cells without counts do not change the pseudobulk",
      identical(symbol_pseudobulk(empty_cell, row_data), summed))
check("genes without counts are dropped",
      identical(names(symbol_pseudobulk(counts[, "b", drop = FALSE], row_data)), "B"))

against <- reference_agreement(c(g1 = 10, g2 = 20, g3 = 30),
    list(list(aligner = "starsolo", max_errors = "0", pseudobulk = c(g1 = 20, g2 = 40, g3 = 60)),
         list(aligner = "starsolo", max_errors = "1", pseudobulk = c(g1 = 20, g2 = 40))))
check("agreement with the reference uses the genes every run holds",
      all(against$n_genes == 2))
check("a doubled yield keeps Spearman and scaled MARD and raises MARD",
      isTRUE(all.equal(against$spearman[1], 1)) && isTRUE(all.equal(against$scaled_mard[1], 0)) &&
      isTRUE(all.equal(against$mard[1], 2 / 3)))
