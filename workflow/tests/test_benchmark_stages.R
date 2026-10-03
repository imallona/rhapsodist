#!/usr/bin/env Rscript
## Self-contained tests for workflow/src/benchmark_stages.R.
## Run with:
##   cd workflow && Rscript tests/test_benchmark_stages.R

script_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
        "."
    }
})
source(file.path(script_dir, "..", "src", "benchmark_stages.R"))
stage_table <- read_stage_table(file.path(script_dir, "..", "data", "benchmark_stages.tsv"))
cwl_stage_table <- read.delim(file.path(script_dir, "..", "data", "sbg_cwl_stages.tsv"),
                              stringsAsFactors = FALSE)

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

stems <- c("star_indexing", "star_s1", "r_sce_generation_s1_star",
           "r_sce_generation_s1_sbg", "alevin_s1_descriptive_report",
           "standardize_cb_umis_s1", "s1_kallisto_bus", "sbg_cwl_s1",
           "s1_sbg_split_sampletags", "mouse_sampletags__index", "made_up_step")
classified <- classify_benchmarks(stems, stage_table)
row <- function(stem) classified[classified$stem == stem, ]

check("one row per stem, in order", identical(classified$stem, stems))
check("a fixed name wins over a sample pattern",
      row("star_indexing")$stage == "reference")
check("a sample pattern matches", row("star_s1")$stage == "alignment and quantification")
check("the aligner suffix sets the pipeline",
      row("r_sce_generation_s1_star")$pipeline == "starsolo" &&
      row("r_sce_generation_s1_sbg")$pipeline == "sbg")
check("an aligner placeholder sets the pipeline",
      row("alevin_s1_descriptive_report")$pipeline == "alevin" &&
      row("s1_sbg_split_sampletags")$pipeline == "sbg")
check("the step keeps its placeholders",
      row("alevin_s1_descriptive_report")$step == "{aligner}_{sample}_descriptive_report")
check("a step used by several aligners is shared",
      row("standardize_cb_umis_s1")$pipeline == "shared" &&
      row("standardize_cb_umis_s1")$aligners == "starsolo,kallisto,alevin")
check("a species placeholder matches", row("mouse_sampletags__index")$stage == "sample tags")
check("an unknown stem is unassigned", row("made_up_step")$stage == "unassigned")

check("named CWL steps get their stage",
      identical(classify_cwl_steps(c("QualCLAlign_RNA", "GetDataTable", "Metrics"), cwl_stage_table),
                c("alignment and quantification", "matrix building", "reporting")))
check("settings, other modalities and unknown CWL steps",
      identical(classify_cwl_steps(c("Bam_Settings", "VDJ_Compile_Results", "QualCLAlign_ATAC", "NewStep"),
                                   cwl_stage_table),
                c("preprocessing", "other modalities", "other modalities", "unassigned")))

seconds <- c(100, 60, 10, 12, 5, 30, 40, 500, 3, 2, 1)
times <- stage_times(classified, seconds, c("starsolo", "kallisto", "sbg"))
value <- function(tab, pipeline, stage) {
    hit <- tab$seconds[tab$pipeline == pipeline & tab$stage == stage]
    if (length(hit) == 0) 0 else hit
}
check("a shared step counts for each aligner using it",
      value(times, "starsolo", "preprocessing") == 30 &&
      value(times, "kallisto", "preprocessing") == 30 &&
      value(times, "sbg", "preprocessing") == 0)
check("an aligner that was not run has no rows", !"alevin" %in% times$pipeline)
check("the CWL rule counts as one stage without step times",
      value(times, "sbg", "alignment and quantification") == 500)

cwl_steps <- data.frame(sample = "s1",
                        step = c("QualCLAlign_RNA", "AnnotateMolecules", "GetDataTable", "Metrics"),
                        seconds = c(300, 100, 50, 20))
split_times <- stage_times(classified, seconds, c("starsolo", "kallisto", "sbg"),
                           cwl_steps, cwl_stage_table)
check("CWL steps replace the CWL rule",
      value(split_times, "sbg", "alignment and quantification") == 400 &&
      value(split_times, "sbg", "matrix building") == 50 + 12 &&
      value(split_times, "sbg", "reporting") == 20)
check("the rest of the CWL rule time is overhead",
      value(split_times, "sbg", "CWL overhead") == 30)
check("other aligners are unchanged by the CWL steps",
      value(split_times, "starsolo", "alignment and quantification") == 60)

cat("\nAll benchmark_stages tests passed.\n")
