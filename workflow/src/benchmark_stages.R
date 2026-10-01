## Assigns benchmark files and BD CWL steps to a pipeline and a stage, from
## workflow/data/benchmark_stages.tsv and workflow/data/sbg_cwl_stages.tsv.

BENCHMARK_ALIGNERS <- c("starsolo", "kallisto", "alevin", "sbg")

## Stages compared between aligners, in plotting order. Download, install,
## reference and simulation steps run once and are left out.
COMPARED_STAGES <- c("preprocessing", "alignment and quantification",
                     "matrix building", "sample tags", "reporting",
                     "extra outputs", "unassigned")

read_stage_table <- function(path) {
    stages <- read.delim(path, stringsAsFactors = FALSE)
    with_aligner <- grepl("{aligner}", stages$pattern, fixed = TRUE)
    expanded <- do.call(rbind, lapply(BENCHMARK_ALIGNERS, function(aligner) {
        rows <- stages[with_aligner, ]
        rows$step <- rows$pattern
        rows$pattern <- gsub("{aligner}", aligner, rows$pattern, fixed = TRUE)
        rows$aligners <- aligner
        rows
    }))
    plain <- stages[!with_aligner, ]
    plain$step <- plain$pattern
    rbind(plain, expanded)
}

pattern_to_regex <- function(pattern) {
    regex <- gsub("{sample}", "(.+)", pattern, fixed = TRUE)
    regex <- gsub("{species}", "(human|mouse)", regex, fixed = TRUE)
    paste0("^", regex, "$")
}

## Characters of a pattern outside its placeholders. The pattern with the most
## of them is taken when several match one file.
literal_length <- function(pattern) nchar(gsub("\\{[a-z]+\\}", "", pattern))

## One row per benchmark file stem: step, pipeline, aligners and stage. A stem
## matching no pattern gets the stage "unassigned".
classify_benchmarks <- function(stems, stage_table) {
    regexes <- vapply(stage_table$pattern, pattern_to_regex, character(1))
    literals <- literal_length(stage_table$pattern)
    rows <- lapply(stems, function(stem) {
        hits <- which(vapply(regexes, grepl, logical(1), x = stem))
        if (length(hits) == 0) {
            return(data.frame(stem = stem, step = stem, pipeline = "shared",
                              aligners = "", stage = "unassigned"))
        }
        best <- stage_table[hits[which.max(literals[hits])], ]
        single <- !grepl(",", best$aligners, fixed = TRUE)
        data.frame(stem = stem, step = best$step,
                   pipeline = if (single) best$aligners else "shared",
                   aligners = best$aligners, stage = best$stage)
    })
    do.call(rbind, rows)
}

## Stage of each BD CWL step; settings steps count as preprocessing and steps
## of the VDJ and ATAC branches as other modalities.
classify_cwl_steps <- function(steps, cwl_stage_table) {
    stage <- cwl_stage_table$stage[match(steps, cwl_stage_table$step)]
    stage[is.na(stage) & grepl("_Settings$|^Start_Time$|^Version$|^GetMachineResources$", steps)] <- "preprocessing"
    stage[is.na(stage) & grepl("VDJ|ATAC|Peak_Annotation", steps)] <- "other modalities"
    stage[is.na(stage)] <- "unassigned"
    stage
}

## Seconds per aligner and stage. A step used by several aligners is counted
## for each of them. When cwl_steps is given (columns sample, step, seconds),
## the single BD CWL rule is replaced by its steps, and the rest of the rule
## time is reported as "CWL overhead".
stage_times <- function(classified, seconds, aligners_run,
                        cwl_steps = NULL, cwl_stage_table = NULL) {
    classified$seconds <- seconds
    per_aligner <- do.call(rbind, lapply(aligners_run, function(aligner) {
        used <- vapply(strsplit(classified$aligners, ",", fixed = TRUE),
                       function(a) aligner %in% a, logical(1))
        rows <- classified[used, c("step", "stage", "seconds")]
        if (nrow(rows) == 0) return(NULL)
        cbind(pipeline = aligner, rows)
    }))
    if (!is.null(cwl_steps) && nrow(cwl_steps) > 0) {
        is_cwl <- per_aligner$pipeline == "sbg" & per_aligner$step == "sbg_cwl_{sample}"
        rule_seconds <- sum(per_aligner$seconds[is_cwl])
        step_rows <- data.frame(
            pipeline = "sbg", step = cwl_steps$step,
            stage = classify_cwl_steps(cwl_steps$step, cwl_stage_table),
            seconds = cwl_steps$seconds)
        overhead <- data.frame(
            pipeline = "sbg", step = "CWL overhead", stage = "CWL overhead",
            seconds = max(0, rule_seconds - sum(cwl_steps$seconds, na.rm = TRUE)))
        per_aligner <- rbind(per_aligner[!is_cwl, ], step_rows, overhead)
    }
    aggregate(seconds ~ pipeline + stage, data = per_aligner, FUN = sum)
}
