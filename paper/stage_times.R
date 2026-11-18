## Time per stage and aligner of one run, from its benchmark files and the
## stage tables under workflow/data, for the paper figures.

stage_times_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "paper"
})
source(file.path(stage_times_dir, "..", "workflow", "src", "benchmark_stages.R"))
STAGE_DATA_DIR <- file.path(stage_times_dir, "..", "workflow", "data")

## Install rules run once; report rules do not run a tool. Both are left out
## of the time and memory of an aligner. Matched by benchmark file name.
is_install_rule <- function(file) grepl("_install\\.txt$", file)
is_report_rule <- function(file) grepl("_report\\.txt$", file)
is_tool_rule <- function(file) !is_install_rule(file) & !is_report_rule(file)

## Seconds of each benchmark file in bench_dir, by stem.
benchmark_seconds <- function(bench_dir) {
    files <- list.files(bench_dir, pattern = "\\.txt$", full.names = TRUE)
    ## the benchmarks report times itself after rendering; an earlier render is left out
    files <- files[basename(files) != "benchmarks_report.txt"]
    seconds <- vapply(files, function(fn) {
        bm <- tryCatch(read.delim(fn, nrows = 1), error = function(e) NULL)
        if (is.null(bm) || !"s" %in% names(bm)) NA_real_ else as.numeric(bm$s[1])
    }, numeric(1))
    data.frame(stem = sub("\\.txt$", "", basename(files)), seconds = seconds,
               stringsAsFactors = FALSE)[!is.na(seconds), ]
}

## BD CWL step times written next to the benchmark files, or NULL.
cwl_step_seconds <- function(bench_dir) {
    files <- list.files(bench_dir, pattern = "^sbg_cwl_steps_.*\\.tsv$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rows <- do.call(rbind, lapply(files, function(fn) {
        steps <- read.delim(fn, stringsAsFactors = FALSE)
        steps <- steps[!is.na(steps$seconds), ]
        if (nrow(steps) == 0) return(NULL)
        data.frame(sample = sub("^sbg_cwl_steps_(.*)\\.tsv$", "\\1", basename(fn)),
                   step = steps$step,
                   seconds = exclusive_seconds(steps$start, steps$end, steps$seconds),
                   stringsAsFactors = FALSE)
    }))
    if (is.null(rows) || nrow(rows) == 0) NULL else rows
}

## One row per aligner and compared stage: seconds and minutes.
stage_time_table <- function(bench_dir, aligners, data_dir = STAGE_DATA_DIR) {
    stage_table <- read_stage_table(file.path(data_dir, "benchmark_stages.tsv"))
    cwl_stage_table <- read.delim(file.path(data_dir, "sbg_cwl_stages.tsv"),
                                  stringsAsFactors = FALSE)
    bm <- benchmark_seconds(bench_dir)
    classified <- classify_benchmarks(bm$stem, stage_table)
    by_stage <- stage_times(classified, bm$seconds, aligners,
                            cwl_step_seconds(bench_dir), cwl_stage_table)
    kept <- setdiff(c(COMPARED_STAGES, "CWL overhead"), "reporting")
    by_stage <- by_stage[by_stage$stage %in% kept & by_stage$seconds > 0, ]
    by_stage$stage <- factor(by_stage$stage, levels = kept)
    by_stage$pipeline <- factor(by_stage$pipeline, levels = aligners)
    by_stage <- by_stage[order(by_stage$pipeline, by_stage$stage), ]
    by_stage$minutes <- by_stage$seconds / 60
    rownames(by_stage) <- NULL
    by_stage
}

## Peak resident memory (GiB) per aligner over the rules of its compared
## stages, as for the time per stage.
aligner_peak_memory <- function(bench_dir, aligners, data_dir = STAGE_DATA_DIR) {
    stage_table <- read_stage_table(file.path(data_dir, "benchmark_stages.tsv"))
    files <- list.files(bench_dir, pattern = "\\.txt$", full.names = TRUE)
    rss <- vapply(files, function(fn) {
        bm <- tryCatch(read.delim(fn, nrows = 1), error = function(e) NULL)
        if (is.null(bm) || !"max_rss" %in% names(bm)) NA_real_ else as.numeric(bm$max_rss[1])
    }, numeric(1))
    classified <- classify_benchmarks(sub("\\.txt$", "", basename(files)), stage_table)
    kept <- classified$stage %in% setdiff(COMPARED_STAGES, "reporting") & !is.na(rss)
    rows <- lapply(aligners, function(aligner) {
        used <- kept & vapply(strsplit(classified$aligners, ",", fixed = TRUE),
                              function(a) aligner %in% a, logical(1))
        if (!any(used)) return(NULL)
        data.frame(pipeline = aligner, peak_rss_gb = max(rss[used]) / 1024)
    })
    do.call(rbind, rows)
}

## The machine record of a run as key and value columns, or NULL when absent.
read_run_info <- function(working_dir) {
    fn <- file.path(working_dir, "run_info.tsv")
    if (!file.exists(fn)) return(NULL)
    read.delim(fn, stringsAsFactors = FALSE, colClasses = "character")
}
