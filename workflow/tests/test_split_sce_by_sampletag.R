#!/usr/bin/env Rscript
## Tests for split_sce_by_sampletag.R and sce_io.R: partition logic, memory/delayed
## backend equivalence, source-count match, h5ad round-trip, and a backend timing table.
## Needs the r_bioc stack. Run: cd workflow && Rscript tests/test_split_sce_by_sampletag.R

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(HDF5Array)
    library(DelayedArray)
    library(data.table)
    library(Matrix)
})

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
split_script <- normalizePath(file.path(script_dir, "..", "src", "split_sce_by_sampletag.R"))
source(file.path(script_dir, "..", "src", "sce_io.R"))

check <- function(label, condition) {
    if (!isTRUE(condition)) stop(sprintf("FAIL: %s", label), call. = FALSE)
    cat(sprintf("ok  - %s\n", label))
}

## Build a synthetic HDF5-backed aligner SCE plus a matching demux table. A fraction of
## cells are non-singlet (status doublet or NA tag) and must be excluded from the splits.
make_case <- function(dir, n_genes, n_cells, n_tags, seed = 1) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    set.seed(seed)
    counts <- as(matrix(rpois(n_genes * n_cells, lambda = 0.5),
                        nrow = n_genes, ncol = n_cells), "CsparseMatrix")
    genes <- sprintf("g%04d", seq_len(n_genes))
    cells <- sprintf("cb%05d", seq_len(n_cells))
    rownames(counts) <- genes
    colnames(counts) <- cells
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                rowData = DataFrame(name = genes))
    sce_dir <- file.path(dir, "source_hdf5")
    sce_h <- saveHDF5SummarizedExperiment(sce, dir = sce_dir, replace = TRUE)
    sce_rds <- file.path(dir, "source_sce.rds")
    base::saveRDS(sce_h, sce_rds)

    tag_names <- sprintf("tag%d", seq_len(n_tags))
    assigned <- tag_names[(seq_len(n_cells) %% n_tags) + 1]
    status <- rep("highqual", n_cells)
    ## make every 7th cell a non-singlet so the split must drop it
    drop <- seq_len(n_cells) %% 7 == 0
    status[drop] <- "multiplet"
    called <- assigned
    called[drop] <- NA
    demux <- data.frame(cb = cells, status = status, called_tag = called,
                        stringsAsFactors = FALSE)
    demux_path <- file.path(dir, "demux.tsv.gz")
    fwrite(demux, demux_path, sep = "\t")

    singlets <- demux[demux$status %in% c("highqual", "called") & !is.na(demux$called_tag), ]
    list(sce_rds = sce_rds, demux = demux_path, demux_df = demux, tags = tag_names,
         counts = counts, singlets = singlets)
}

run_split <- function(case, out_base, backend, output_format = "sce") {
    args <- c("-q", "--no-save", "--no-restore", "--slave", "-f", split_script, "--args",
              "--sce", case$sce_rds, "--demux", case$demux, "--out_base", out_base,
              "--backend", backend, "--output_format", output_format,
              "--tags", paste(case$tags, collapse = ","))
    t <- system.time(status <- system2("R", args, stdout = NULL, stderr = NULL))
    if (status != 0) stop(sprintf("split failed (backend=%s, format=%s)", backend, output_format))
    t[["elapsed"]]
}

load_split_counts <- function(out_base, label) {
    as.matrix(counts(loadHDF5SummarizedExperiment(file.path(out_base, label))))
}

tmp <- tempfile("split_test_")
dir.create(tmp)

## unit test: split_sce_singlets partition logic
case <- make_case(file.path(tmp, "unit"), n_genes = 30, n_cells = 60, n_tags = 4)
sce_mem <- realize_assays_in_memory(readRDS(case$sce_rds))
splits <- split_sce_singlets(sce_mem, case$singlets, case$tags, case$tags)

check("one split per tag", length(splits) == length(case$tags))
all_cells <- unlist(lapply(splits, colnames), use.names = FALSE)
check("splits are disjoint", !any(duplicated(all_cells)))
check("union of splits equals the singlet cells present in the SCE",
      setequal(all_cells, intersect(case$singlets$cb, colnames(sce_mem))))
check("non-singlet cells are excluded",
      length(intersect(all_cells, case$demux_df$cb[is.na(case$demux_df$called_tag)])) == 0)
for (i in seq_along(case$tags)) {
    want <- case$singlets$cb[case$singlets$called_tag == case$tags[i]]
    check(sprintf("tag %s holds exactly its singlet cells", case$tags[i]),
          setequal(colnames(splits[[i]]), intersect(want, colnames(sce_mem))))
}

## integration: memory and delayed backends agree, and match the source counts
mem_base <- file.path(tmp, "mem")
del_base <- file.path(tmp, "del")
invisible(run_split(case, mem_base, "memory"))
invisible(run_split(case, del_base, "delayed"))

for (tag in case$tags) {
    cm <- load_split_counts(mem_base, tag)
    cd <- load_split_counts(del_base, tag)
    check(sprintf("memory and delayed backends give identical cells for %s", tag),
          identical(colnames(cm), colnames(cd)))
    check(sprintf("memory and delayed backends give identical counts for %s", tag),
          identical(cm, cd))
    expected <- as.matrix(case$counts[, colnames(cm), drop = FALSE])
    check(sprintf("split counts match the source for %s", tag),
          identical(cm, expected))
}

## h5ad round-trip (only when anndataR is installed)
if (requireNamespace("anndataR", quietly = TRUE)) {
    both_base <- file.path(tmp, "both")
    run_split(case, both_base, "memory", output_format = "both")
    tag <- case$tags[1]
    h5ad_path <- file.path(both_base, tag, "adata.h5ad")
    check("h5ad split file is written under output_format both", file.exists(h5ad_path))
    back <- anndataR::read_h5ad(h5ad_path, as = "SingleCellExperiment")
    sce_counts <- load_split_counts(both_base, tag)
    back_counts <- as.matrix(counts(back))
    ## anndataR keeps genes x cells via the SCE reader; align dimnames before comparing
    back_counts <- back_counts[rownames(sce_counts), colnames(sce_counts), drop = FALSE]
    check("h5ad counts round-trip to the SCE split counts",
          identical(back_counts, sce_counts))
} else {
    cat("skip - anndataR not installed; h5ad round-trip not checked\n")
}

## A DelayedArray seed that counts how often the source is read (extract_array). The
## delayed backend reads it once per tag, the memory backend once in total, so the read
## count is asserted below. utime is printed but not asserted: this seed is in-memory, so
## it does not show the on-disk HDF5 cost that makes the delayed backend slow.
setClass("CountingSeed", representation(data = "matrix", counter = "environment"))
setMethod("dim", "CountingSeed", function(x) dim(x@data))
setMethod("dimnames", "CountingSeed", function(x) dimnames(x@data))
setMethod("extract_array", "CountingSeed", function(x, index) {
    x@counter$reads <- x@counter$reads + 1L
    DelayedArray::extract_array(x@data, index)
})

split_counting_source <- function(n_genes, n_cells, n_tags, backend, out_base) {
    set.seed(7)
    m <- matrix(rpois(n_genes * n_cells, 0.5), n_genes, n_cells)
    rownames(m) <- sprintf("g%04d", seq_len(n_genes))
    colnames(m) <- sprintf("cb%05d", seq_len(n_cells))
    counter <- new.env(); counter$reads <- 0L
    sce <- SingleCellExperiment(assays = list(
        counts = DelayedArray(new("CountingSeed", data = m, counter = counter))))
    tags <- sprintf("tag%d", seq_len(n_tags))
    singlets <- data.frame(cb = colnames(m),
                           called_tag = tags[(seq_len(n_cells) %% n_tags) + 1],
                           stringsAsFactors = FALSE)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
    ut <- system.time({
        if (backend == "memory") sce <- realize_assays_in_memory(sce)
        sp <- split_sce_singlets(sce, singlets, tags, tags)
        for (i in seq_along(tags)) {
            saveHDF5SummarizedExperiment(sp[[i]], dir = file.path(out_base, tags[i]),
                                         replace = TRUE)
        }
    })
    list(reads = counter$reads, utime = ut[["user.self"]] + ut[["sys.self"]])
}

few <- 3; many <- 9
runs <- list(
    mem_few  = split_counting_source(300, 1500, few,  "memory",  file.path(tmp, "rd_mem_few")),
    mem_many = split_counting_source(300, 1500, many, "memory",  file.path(tmp, "rd_mem_many")),
    del_few  = split_counting_source(300, 1500, few,  "delayed", file.path(tmp, "rd_del_few")),
    del_many = split_counting_source(300, 1500, many, "delayed", file.path(tmp, "rd_del_many")))

cat("\nsource pulls (extract_array calls) and CPU time (utime, s):\n")
for (nm in names(runs)) {
    cat(sprintf("  %-9s reads=%2d  utime=%.2f\n", nm, runs[[nm]]$reads, runs[[nm]]$utime))
}

check("memory backend pulls the source a fixed number of times regardless of tag count",
      runs$mem_few$reads == runs$mem_many$reads)
check("delayed backend pulls the source more times as the tag count grows",
      runs$del_many$reads > runs$del_few$reads)
check("delayed pulls scale with the tag count (at least one extra pass per added tag)",
      runs$del_many$reads - runs$del_few$reads >= many - few)
check("memory backend pulls the source far fewer times than delayed at many tags",
      runs$mem_many$reads < runs$del_many$reads)

unlink(tmp, recursive = TRUE)
cat("\nAll sampletag split tests passed.\n")
