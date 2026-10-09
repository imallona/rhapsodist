#!/usr/bin/env Rscript
## Tests for paper/compare_sampletag_calls.R. Run from the repository root:
##   Rscript paper/test_compare_sampletag_calls.R

sys.source("paper/compare_sampletag_calls.R", envir = environment())

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

whitelists <- list(c("A1", "A2", "A3"), c("B1", "B2", "B3"), c("C1", "C2", "C3"))
check("index 1 is the first entry of each segment",
      cell_index_to_barcode(1, whitelists) == "A1B1C1")
check("the last segment changes fastest",
      identical(cell_index_to_barcode(c(2, 3, 4), whitelists), c("A1B1C2", "A1B1C3", "A1B2C1")))
check("the last index is the last entry of each segment",
      cell_index_to_barcode(27, whitelists) == "A3B3C3")

deposited <- data.table(deposited_tag = c("t1", "t1", "t1", "t2", "t2"),
                        cb = c("a", "b", "c", "d", "e"))
demux <- data.table(cb = c("a", "b", "c", "d"),
                    called_tag = c("t1", "t2", NA, NA),
                    status = c("highqual", "called", "multiplet", "undetermined"))
calls <- classify_calls(deposited, demux)
outcome <- setNames(calls$outcome, calls$cb)
check("every outcome is recognised",
      identical(unname(outcome[c("a", "b", "c", "d", "e")]),
                c("same tag", "other tag", "multiplet", "undetermined", "absent")))

cat("\nAll compare_sampletag_calls tests passed.\n")
