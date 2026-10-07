## Time per stage and aligner of one run, from its benchmark files and the
## stage tables under workflow/data, for the paper figures.

stage_times_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "paper"
})
source(file.path(stage_times_dir, "..", "workflow", "src", "benchmark_stages.R"))
STAGE_DATA_DIR <- file.path(stage_times_dir, "..", "workflow", "data")

## Seconds of each benchmark file in bench_dir, by stem.
benchmark_seconds <- function(bench_dir) {
    files <- list.files(bench_dir, pattern = "\\.txt$", full.names = TRUE)
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
        data.frame(sample = sub("^sbg_cwl_steps_(.*)\\.tsv$", "\\1", basename(fn)),
                   step = steps$step, seconds = steps$seconds, stringsAsFactors = FALSE)
    }))
    if (nrow(rows) == 0) NULL else rows
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
    kept <- c(COMPARED_STAGES, "CWL overhead")
    by_stage <- by_stage[by_stage$stage %in% kept & by_stage$seconds > 0, ]
    by_stage$stage <- factor(by_stage$stage, levels = kept)
    by_stage$pipeline <- factor(by_stage$pipeline, levels = aligners)
    by_stage <- by_stage[order(by_stage$pipeline, by_stage$stage), ]
    by_stage$minutes <- by_stage$seconds / 60
    rownames(by_stage) <- NULL
    by_stage
}

## The machine record of a run as key and value columns, or NULL when absent.
read_run_info <- function(working_dir) {
    fn <- file.path(working_dir, "run_info.tsv")
    if (!file.exists(fn)) return(NULL)
    read.delim(fn, stringsAsFactors = FALSE, colClasses = "character")
}
