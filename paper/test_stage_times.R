#!/usr/bin/env Rscript
## Tests for paper/stage_times.R. Run from the repository root:
##   Rscript paper/test_stage_times.R

source("paper/stage_times.R")

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

write_benchmark <- function(dir, stem, seconds) {
    writeLines(c("s\th:m:s\tmax_rss", sprintf("%s\t0:00:00\t1", seconds)),
               file.path(dir, paste0(stem, ".txt")))
}

run <- tempfile("run")
bench <- file.path(run, "benchmarks")
dir.create(bench, recursive = TRUE)
write_benchmark(bench, "standardize_cb_umis_s1", 60)
write_benchmark(bench, "star_s1", 600)
write_benchmark(bench, "s1_bustools_count", 120)
write_benchmark(bench, "sbg_cwl_s1", 1000)
write_benchmark(bench, "download_genome", 500)
writeLines(c("step\tstart\tend\tseconds\tstatus",
             "QualCLAlign_RNA\ta\tb\t700\tsuccess",
             "Metrics\ta\tb\tNA\tfailed"),
           file.path(bench, "sbg_cwl_steps_s1.tsv"))

aligners <- c("starsolo", "kallisto", "alevin", "sbg")
tab <- stage_time_table(bench, aligners)

check("one row per aligner and stage with time",
      !any(duplicated(tab[, c("pipeline", "stage")])))
check("a shared step is counted for each aligner that uses it",
      all(tab$seconds[tab$stage == "preprocessing" & tab$pipeline %in% c("starsolo", "kallisto")] == 60))
check("alignment time per aligner",
      tab$seconds[tab$pipeline == "starsolo" & tab$stage == "alignment and quantification"] == 600 &&
      tab$seconds[tab$pipeline == "kallisto" & tab$stage == "alignment and quantification"] == 120)
check("the CWL rule is split into its steps and an overhead",
      tab$seconds[tab$pipeline == "sbg" & tab$stage == "alignment and quantification"] == 700 &&
      tab$seconds[tab$pipeline == "sbg" & tab$stage == "CWL overhead"] == 300)
check("download steps are left out", !"download" %in% as.character(tab$stage))
check("an aligner with only shared steps has only those stages",
      identical(as.character(tab$stage[tab$pipeline == "alevin"]), "preprocessing"))
check("minutes derive from seconds", all(abs(tab$minutes * 60 - tab$seconds) < 1e-9))
check("aligners are ordered as given",
      identical(levels(tab$pipeline), aligners))

empty <- file.path(run, "empty_benchmarks")
dir.create(empty)
writeLines("step\tstart\tend\tseconds\tstatus", file.path(empty, "sbg_cwl_steps_s2.tsv"))
write_benchmark(empty, "sbg_cwl_s2", 50)
check("a header-only CWL step file gives no steps", is.null(cwl_step_seconds(empty)))
only_rule <- stage_time_table(empty, "sbg")
check("the CWL rule time is kept when no step was timed",
      only_rule$seconds[only_rule$stage == "alignment and quantification"] == 50)

check("missing run info is NULL", is.null(read_run_info(run)))
writeLines(c("key\tvalue", "cpu_model\tAMD EPYC 7763", "n_cpus\t128"),
           file.path(run, "run_info.tsv"))
info <- read_run_info(run)
check("run info is read as text", identical(info$value[info$key == "n_cpus"], "128"))

cat("all stage time tests passed\n")
