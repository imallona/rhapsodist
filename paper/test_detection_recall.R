#!/usr/bin/env Rscript
## Tests for paper/detection_recall.R. Run from the repository root:
##   Rscript paper/test_detection_recall.R

sys.source("paper/detection_recall.R", envir = environment())

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

truth <- matrix(c(0, 1, 1, 4,
                  0, 0, 2, 12),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("g1", "g2"), c("a", "b", "c", "d")))
observed <- matrix(c(1, 1, 0, 3, 9,
                     0, 0, 2, 0, 9,
                     5, 5, 5, 5, 9),
                   nrow = 3, byrow = TRUE,
                   dimnames = list(c("g1", "g2", "g3"), c("a", "b", "c", "d", "e")))
result <- detection_by_true_count(truth, observed)
fraction <- setNames(result$detected_fraction, result$true_count)
entries <- setNames(result$n_entries, result$true_count)

check("genes and cells absent from the truth are ignored", sum(entries) == 8)
check("true zeros give the false detection fraction", fraction[["0"]] == 1 / 3)
check("recall at a true count of 1", fraction[["1"]] == 1 / 2)
check("recall at a true count of 2", fraction[["2"]] == 1)
check("bins group the higher counts", entries[["3-5"]] == 1 && entries[[">10"]] == 1)
check("an undetected high count has recall 0", fraction[[">10"]] == 0)
check("an empty bin has no fraction", is.na(fraction[["6-10"]]) && entries[["6-10"]] == 0)
