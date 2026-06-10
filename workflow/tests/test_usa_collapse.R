#!/usr/bin/env Rscript
## Tests for the USA assay collapse in workflow/src/generate_sce_alevin.R.
## alevin-fry names gene columns by splicing status: spliced is the bare gene_id,
## unspliced is <gene_id>-U, ambiguous is <gene_id>-A (no -S suffix). The collapse
## must route each status to its own assay and set counts = spliced + ambiguous.
## Run with:
##   cd workflow && Rscript tests/test_usa_collapse.R

suppressPackageStartupMessages(library(Matrix))

script_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    f_flag <- which(ca == "-f")
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else if (length(f_flag) > 0 && length(ca) > f_flag[1]) {
        dirname(normalizePath(ca[f_flag[1] + 1L]))
    } else {
        "."
    }
})

## evaluate only the collapse_usa_assays definition from the real script, so the
## test exercises the committed function rather than a copy that can drift
src_path <- file.path(script_dir, "..", "src", "generate_sce_alevin.R")
exprs <- parse(src_path)
for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(e[[2]], as.name("collapse_usa_assays"))) {
        eval(e)
    }
}

check <- function(label, condition) {
    if (!isTRUE(condition)) stop(sprintf("FAIL: %s", label), call. = FALSE)
    cat(sprintf("ok  - %s\n", label))
}

check("collapse_usa_assays was extracted from the script", exists("collapse_usa_assays"))

## gene1 has all three statuses; gene2 is spliced-only (no -U/-A rows), which
## happens for genes with no intronic or ambiguous reads
rn <- c("gene1", "gene1-U", "gene1-A", "gene2")
m <- Matrix(0, nrow = length(rn), ncol = 2, sparse = TRUE,
            dimnames = list(rn, c("c1", "c2")))
m["gene1", "c1"] <- 5; m["gene1-U", "c1"] <- 2; m["gene1-A", "c1"] <- 1
m["gene2", "c2"] <- 7

a <- collapse_usa_assays(m)

check("returns counts, spliced, unspliced, ambiguous",
      identical(names(a), c("counts", "spliced", "unspliced", "ambiguous")))
check("collapses to one row per base gene id",
      identical(rownames(a$counts), c("gene1", "gene2")))
check("spliced takes the bare gene_id column",
      a$spliced["gene1", "c1"] == 5)
check("unspliced takes the -U column",
      a$unspliced["gene1", "c1"] == 2)
check("ambiguous takes the -A column",
      a$ambiguous["gene1", "c1"] == 1)
check("counts equals spliced plus ambiguous",
      all(a$counts == a$spliced + a$ambiguous))
check("counts for a mixed gene is spliced+ambiguous (5+1)",
      a$counts["gene1", "c1"] == 6)
check("spliced-only gene keeps its spliced count and zero U/A",
      a$counts["gene2", "c2"] == 7 &&
      a$unspliced["gene2", "c2"] == 0 &&
      a$ambiguous["gene2", "c2"] == 0)
check("no signal is lost: sum of S+U+A equals raw input sum",
      sum(a$spliced) + sum(a$unspliced) + sum(a$ambiguous) == sum(m))

cat("\nAll USA collapse tests passed.\n")
