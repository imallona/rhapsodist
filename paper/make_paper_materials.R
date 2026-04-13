#!/usr/bin/env Rscript
## make_paper_materials.R
##
## Reads the RDS outputs produced by the main rhapsodist pipeline and writes
## manuscript-ready CSV and PDF artefacts under <working_dir>/paper/<sample>/
## (or <working_dir>/paper/ for benchmark-only outputs).
##
## Usage:
##   Rscript paper/make_paper_materials.R --kind simulation \
##       --working_dir output/simul --sample simulated \
##       --aligners starsolo,kallisto,alevin,sbg --n_expected_cells 1000
##
##   Rscript paper/make_paper_materials.R --kind biology \
##       --working_dir output/sendoel2024 --sample sample_16_wta_p60 \
##       --aligners starsolo,kallisto,alevin
##
##   Rscript paper/make_paper_materials.R --kind benchmarks \
##       --working_dir output/simul --n_expected_cells 1000
##
## Inputs are the RDS files already saved by 02_comparison.Rmd,
## 04_biology.Rmd, and the benchmark .txt files written by Snakemake.
## Nothing here depends on the production reports.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(ggrepel)
    library(Matrix)
    library(SingleCellExperiment)
})

`%||%` = function(a, b) if (!is.null(a)) a else b

script_dir = tryCatch({
    args_file = commandArgs(trailingOnly = FALSE)
    f = sub("^--file=", "", grep("^--file=", args_file, value = TRUE))
    if (length(f) == 0) getwd() else dirname(normalizePath(f))
}, error = function(e) getwd())
source(file.path(script_dir, "export_helpers.R"))

parse_args = function() {
    p = ArgumentParser()
    p$add_argument("--kind", required = TRUE,
                   choices = c("simulation", "biology", "benchmarks"))
    p$add_argument("--working_dir", required = TRUE)
    p$add_argument("--sample", default = "")
    p$add_argument("--aligners", default = "starsolo,kallisto,alevin")
    p$add_argument("--has_sbg", default = "false")
    p$add_argument("--n_expected_cells", type = "integer", default = 0L)
    p$add_argument("--markers_file", default = "")
    p$parse_args()
}

aligner_colours = c(starsolo = "#E69F00", kallisto = "#56B4E9",
                    alevin = "#009E73", sbg = "#CC79A7",
                    truth = "#000000")

## Simulation kind. Reads per-pipeline SCEs plus the simulation truth (barcodes
## and matrix market format) and produces Figs 1 to 3 plus S1, S2, Table 1.
run_simulation = function(opt) {
    wd = opt$working_dir
    samp = opt$sample
    aligners = strsplit(opt$aligners, ",")[[1]]
    if (isTRUE(as.logical(opt$has_sbg))) aligners = c(aligners, "sbg")
    pdir = setup_paper_dir(wd, samp)

    sces = lapply(setNames(aligners, aligners), function(pipe) {
        fn = file.path(wd, pipe, samp, paste0(samp, "_", pipe, "_sce.rds"))
        if (!file.exists(fn)) {
            message("SCE missing: ", fn); return(NULL)
        }
        readRDS(fn)
    })
    sces = Filter(Negate(is.null), sces)

    ## Cell counts.
    counts_df = data.table(
        pipeline = names(sces),
        n_cells = vapply(sces, ncol, integer(1)))
    pp_save_csv(counts_df, pdir, "sim_cell_counts")
    p_counts = ggplot(counts_df, aes(pipeline, n_cells, fill = pipeline)) +
        geom_col() + geom_text(aes(label = n_cells), vjust = -0.3, size = 3) +
        scale_fill_manual(values = aligner_colours) +
        theme_bw() + theme(legend.position = "none") +
        labs(x = NULL, y = "cells recovered")
    pp_save_pdf(p_counts, pdir, "sim_cell_counts", width = 4, height = 3)

    ## Ground-truth recall and precision if truth barcodes are present.
    bc_file = file.path(wd, "simulate", "cell_barcodes.txt")
    if (file.exists(bc_file)) {
        truth_bc = readLines(bc_file)
        gt = rbindlist(lapply(names(sces), function(pipe) {
            cb = colnames(sces[[pipe]])
            data.table(pipeline = pipe,
                       n_recovered = length(cb),
                       n_truth = length(truth_bc),
                       recall = mean(cb %in% truth_bc),
                       precision = sum(cb %in% truth_bc) / length(truth_bc))
        }))
        pp_save_csv(gt, pdir, "sim_ground_truth")

        ## Upset of cell barcodes including truth.
        bc_lists = c(lapply(sces, colnames), list(truth = truth_bc))
        if (requireNamespace("UpSetR", quietly = TRUE) && length(bc_lists) >= 2) {
            pp_save_base_pdf(
                quote(print(UpSetR::upset(UpSetR::fromList(bc_lists),
                                          order.by = "freq", nsets = length(bc_lists),
                                          text.scale = 1.2))),
                pdir, "sim_bc_upset", width = 6, height = 4)
            set_sizes = data.table(set = names(bc_lists),
                                   size = vapply(bc_lists, length, integer(1)))
            pp_save_csv(set_sizes, pdir, "sim_bc_upset_sizes")
        }
    }

    ## Pseudobulk MARD including truth, plus bootstrap CIs.
    truth_mex = file.path(wd, "simulate", "true_mex")
    pb_list = lapply(sces, function(sce) rowSums(counts(sce)))

    truth_pb = NULL
    if (dir.exists(truth_mex)) {
        mat_fn = list.files(truth_mex, pattern = "matrix\\.mtx(\\.gz)?$", full.names = TRUE)[1]
        feat_fn = list.files(truth_mex, pattern = "features\\.tsv(\\.gz)?$", full.names = TRUE)[1]
        if (!is.na(mat_fn) && !is.na(feat_fn)) {
            m = readMM(mat_fn)
            feat = read.table(feat_fn, sep = "\t", stringsAsFactors = FALSE)
            rownames(m) = feat$V1
            truth_pb = rowSums(m)
            pb_list$truth = truth_pb
        }
    }

    shared = Reduce(intersect, lapply(pb_list, names))
    if (length(shared) > 0 && length(pb_list) >= 2) {
        pb_mat = do.call(cbind, lapply(pb_list, function(v) v[shared]))
        mard_stat = function(x, y) {
            d = (x + y) / 2
            mean(abs(x - y)[d > 0] / d[d > 0], na.rm = TRUE)
        }
        pairs = expand.grid(a = colnames(pb_mat), b = colnames(pb_mat),
                            stringsAsFactors = FALSE)
        mard_dt = rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
            a = pairs$a[i]; b = pairs$b[i]
            x = pb_mat[, a]; y = pb_mat[, b]
            m = mard_stat(x, y)
            ci = if (a == b) c(NA, NA) else pp_bootstrap_pair(x, y, mard_stat, n_boot = 1000)
            data.table(pipeline1 = a, pipeline2 = b, mard = m,
                       ci_lo = ci[1], ci_hi = ci[2])
        }))
        pp_save_csv(mard_dt, pdir, "sim_pseudobulk_mard")

        p_mard = ggplot(mard_dt, aes(pipeline1, pipeline2, fill = mard * 100)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.1f%%", mard * 100)), size = 3) +
            scale_fill_distiller(palette = "RdYlGn", direction = -1, limits = c(0, NA)) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
            labs(x = NULL, y = NULL, fill = "MARD (%)")
        pp_save_pdf(p_mard, pdir, "sim_pseudobulk_mard", width = 5, height = 4)

        ## Per-aligner sanity table versus truth.
        if (!is.null(truth_pb)) {
            vs_truth = mard_dt[pipeline2 == "truth" & pipeline1 != "truth"]
            sanity = rbindlist(lapply(names(sces), function(pipe) {
                sce = sces[[pipe]]
                cb = colnames(sce)
                bc_in_truth = if (exists("truth_bc", inherits = FALSE)) {
                    mean(cb %in% truth_bc)
                } else NA_real_
                data.table(
                    pipeline = pipe,
                    n_cells = ncol(sce),
                    n_genes_detected = sum(rowSums(counts(sce)) > 0),
                    median_umi = median(colSums(counts(sce))),
                    mard_vs_truth = vs_truth[pipeline1 == pipe, mard][1],
                    mard_vs_truth_ci_lo = vs_truth[pipeline1 == pipe, ci_lo][1],
                    mard_vs_truth_ci_hi = vs_truth[pipeline1 == pipe, ci_hi][1])
            }))
            pp_save_csv(sanity, pdir, "sim_per_aligner_sanity")
            pp_save_latex_table(sanity, pdir, "sim_per_aligner_sanity",
                caption = "Per-aligner sanity metrics on the simulated dataset.",
                label = "tab:sim_sanity", digits = 3)
        }
    }

    ## UMI distributions per pipeline.
    umi_dt = rbindlist(lapply(names(sces), function(pipe) {
        data.table(pipeline = pipe, umi = colSums(counts(sces[[pipe]])))
    }))
    pp_save_csv(umi_dt[, .(pipeline, median_umi = median(umi),
                           mean_umi = mean(umi)), by = pipeline],
                pdir, "sim_umi_summary")
    p_umi = ggplot(umi_dt, aes(log10(umi + 1), colour = pipeline)) +
        geom_density(linewidth = 0.6) +
        scale_colour_manual(values = aligner_colours) +
        theme_bw() +
        labs(x = "log10(UMI + 1)", y = "density", colour = NULL)
    pp_save_pdf(p_umi, pdir, "sim_umi_distribution", width = 5, height = 3)

    ## Accuracy vs speed if benchmark files exist.
    bench_dir = file.path(wd, "benchmarks")
    if (dir.exists(bench_dir) && !is.null(truth_pb) && exists("pb_mat")) {
        pipe_time = load_pipeline_time(bench_dir, aligners)
        vs_truth = mard_dt[pipeline2 == "truth" & pipeline1 != "truth",
                           .(pipeline = pipeline1, mard, ci_lo, ci_hi)]
        speed_dt = merge(vs_truth, pipe_time, by = "pipeline", all.x = TRUE)
        pp_save_csv(speed_dt, pdir, "sim_accuracy_vs_speed")
        p_acc = ggplot(speed_dt, aes(total_min, mard * 100, colour = pipeline,
                                     label = pipeline)) +
            geom_point(size = 4) +
            geom_errorbar(aes(ymin = ci_lo * 100, ymax = ci_hi * 100),
                          width = 0, linewidth = 0.4) +
            geom_text_repel(size = 3.5, show.legend = FALSE) +
            scale_colour_manual(values = aligner_colours) +
            theme_bw() +
            labs(x = "total pipeline time (min)",
                 y = "MARD vs truth (%)", colour = NULL)
        pp_save_pdf(p_acc, pdir, "sim_accuracy_vs_speed", width = 5, height = 4)
    }

    message("simulation materials written to ", pdir)
}

## Biology kind. Consumes derived biology_*.rds files written by the
## production report and Seurat caches for UMAP plots.
run_biology = function(opt) {
    wd = opt$working_dir
    samp = opt$sample
    aligners = strsplit(opt$aligners, ",")[[1]]
    if (isTRUE(as.logical(opt$has_sbg))) aligners = c(aligners, "sbg")
    pdir = setup_paper_dir(wd, samp)

    biords = function(tag) file.path(wd, paste0(samp, "_biology_", tag, ".rds"))

    ## QC summary (S1 Table).
    if (file.exists(biords("qc"))) {
        qc_df = readRDS(biords("qc"))
        pp_save_csv(qc_df, pdir, "bio_qc")
        pp_save_latex_table(qc_df, pdir, "bio_qc",
            caption = "Per-pipeline QC summary on the real dataset.",
            label = "tab:bio_qc")
        qc_long = melt(as.data.table(qc_df)[, .(pipeline, n_cells,
                                                median_umi, median_genes)],
                       id.vars = "pipeline", variable.name = "metric")
        p_qc = ggplot(qc_long, aes(pipeline, value, fill = pipeline)) +
            geom_col() + facet_wrap(~metric, scales = "free_y") +
            scale_fill_manual(values = aligner_colours) +
            theme_bw() + theme(legend.position = "none") +
            labs(x = NULL, y = NULL)
        pp_save_pdf(p_qc, pdir, "bio_qc_bars", width = 6, height = 3)
    }

    ## Pseudobulk MARD with bootstrap CIs.
    if (file.exists(biords("pseudobulk_mard"))) {
        pbobj = readRDS(biords("pseudobulk_mard"))
        pb_mat = pbobj$pb
        mard_stat = function(x, y) {
            d = (x + y) / 2
            mean(abs(x - y)[d > 0] / d[d > 0], na.rm = TRUE)
        }
        pipes = colnames(pb_mat)
        pairs = expand.grid(a = pipes, b = pipes, stringsAsFactors = FALSE)
        mard_dt = rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
            a = pairs$a[i]; b = pairs$b[i]
            m = mard_stat(pb_mat[, a], pb_mat[, b])
            ci = if (a == b) c(NA, NA) else pp_bootstrap_pair(
                pb_mat[, a], pb_mat[, b], mard_stat, n_boot = 1000)
            data.table(pipeline1 = a, pipeline2 = b, mard = m,
                       ci_lo = ci[1], ci_hi = ci[2])
        }))
        pp_save_csv(mard_dt, pdir, "bio_pseudobulk_mard")
        p = ggplot(mard_dt, aes(pipeline1, pipeline2, fill = mard * 100)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.1f%%", mard * 100)), size = 3) +
            scale_fill_distiller(palette = "RdYlGn", direction = -1, limits = c(0, NA)) +
            theme_bw() +
            labs(x = NULL, y = NULL, fill = "MARD (%)")
        pp_save_pdf(p, pdir, "bio_pseudobulk_mard", width = 4.5, height = 3.5)
    }

    ## Bland-Altman UMI plots. cb_umi holds per-barcode UMIs per pipeline.
    if (file.exists(biords("cb_umi"))) {
        cb_umi = as.data.table(readRDS(biords("cb_umi")))
        pipes = unique(cb_umi$pipeline)
        if (length(pipes) >= 2) {
            pairs = combn(pipes, 2, simplify = FALSE)
            summary_rows = list()
            ba_plots = lapply(pairs, function(pair) {
                d1 = cb_umi[pipeline == pair[1], .(barcode, u1 = total_umi)]
                d2 = cb_umi[pipeline == pair[2], .(barcode, u2 = total_umi)]
                shared = merge(d1, d2, by = "barcode")
                if (nrow(shared) < 10) return(NULL)
                shared[, mean_umi := (log10(u1 + 1) + log10(u2 + 1)) / 2]
                shared[, diff_umi := log10(u1 + 1) - log10(u2 + 1)]
                bias = mean(shared$diff_umi)
                loa = sd(shared$diff_umi) * 1.96
                summary_rows[[paste(pair, collapse = "_vs_")]] <<- data.table(
                    pipeline1 = pair[1], pipeline2 = pair[2],
                    n_shared = nrow(shared),
                    bias = bias, loa_lower = bias - loa, loa_upper = bias + loa)
                ggplot(shared, aes(mean_umi, diff_umi)) +
                    pp_rasterise(geom_point(alpha = 0.1, size = 0.3)) +
                    geom_hline(yintercept = bias, colour = "#D55E00") +
                    geom_hline(yintercept = c(bias - loa, bias + loa),
                               linetype = "dashed", colour = "#D55E00") +
                    theme_bw() +
                    labs(x = "mean log10(UMI + 1)",
                         y = "diff log10(UMI + 1)",
                         title = paste(pair[1], "vs", pair[2]),
                         subtitle = sprintf("n = %d, bias = %.3f",
                                            nrow(shared), bias))
            })
            ba_plots = Filter(Negate(is.null), ba_plots)
            if (length(ba_plots) > 0) {
                pp_save_csv(rbindlist(summary_rows), pdir, "bio_bland_altman")
                pp_save_pdf(wrap_plots(ba_plots, nrow = 1), pdir,
                            "bio_bland_altman",
                            width = 3.2 * length(ba_plots), height = 3.2)
            }
        }
    }

    ## Cluster ARI with bootstrap CIs.
    add_ari_ci = function(cluster_dt, key_a, key_b) {
        pipes = unique(cluster_dt$pipeline)
        if (length(pipes) < 2) return(NULL)
        pairs = combn(pipes, 2, simplify = FALSE)
        rbindlist(lapply(pairs, function(pair) {
            d1 = cluster_dt[pipeline == pair[1],
                            .(barcode, a = get(key_a))]
            d2 = cluster_dt[pipeline == pair[2],
                            .(barcode, b = get(key_b))]
            shared = merge(d1, d2, by = "barcode")
            if (nrow(shared) < 10) return(NULL)
            ari = mclust::adjustedRandIndex(shared$a, shared$b)
            ci = pp_bootstrap_pair(shared$a, shared$b,
                function(x, y) mclust::adjustedRandIndex(x, y),
                n_boot = 1000)
            data.table(pipeline1 = pair[1], pipeline2 = pair[2],
                       n_shared = nrow(shared),
                       ari = ari, ci_lo = ci[1], ci_hi = ci[2])
        }))
    }

    if (file.exists(biords("clusters")) && requireNamespace("mclust", quietly = TRUE)) {
        clusters_dt = as.data.table(readRDS(biords("clusters")))
        cluster_ari = add_ari_ci(clusters_dt, "cluster_prefixed", "cluster_prefixed")
        if (!is.null(cluster_ari)) pp_save_csv(cluster_ari, pdir, "bio_cluster_ari")
        ct_ari = add_ari_ci(clusters_dt, "celltype", "celltype")
        if (!is.null(ct_ari)) pp_save_csv(ct_ari, pdir, "bio_celltype_ari")

        ari_bar = function(dt, name, ylab) {
            if (is.null(dt) || nrow(dt) == 0) return(invisible(NULL))
            dt[, pair := paste(pipeline1, pipeline2, sep = " vs ")]
            p = ggplot(dt, aes(pair, ari)) +
                geom_col(fill = "#0072B2") +
                geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                              width = 0.2, linewidth = 0.4) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
                labs(x = NULL, y = ylab)
            pp_save_pdf(p, pdir, name, width = 4.5, height = 3.2)
        }
        ari_bar(cluster_ari, "bio_cluster_ari", "cluster ARI")
        ari_bar(ct_ari,      "bio_celltype_ari", "cell-type ARI")
    }

    ## UMAP coloured by marker cell type, one panel per pipeline. Requires
    ## Seurat caches. Rasterise the point layer to keep PDF small.
    seu_fns = setNames(file.path(wd, sprintf("%s_biology_%s_seurat.rds", samp,
                                             aligners)), aligners)
    if (all(file.exists(seu_fns)) &&
        requireNamespace("Seurat", quietly = TRUE)) {
        cbbPalette = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        panels = lapply(names(seu_fns), function(pipe) {
            so = readRDS(seu_fns[[pipe]])
            ct = so$marker_celltype %||% rep("none", ncol(so))
            emb = as.data.frame(Seurat::Embeddings(so, "umap"))
            colnames(emb) = c("UMAP_1", "UMAP_2")
            emb$celltype = ct
            ggplot(emb, aes(UMAP_1, UMAP_2, colour = celltype)) +
                pp_rasterise(geom_point(size = 0.3, alpha = 0.8)) +
                scale_colour_manual(values = setNames(
                    rep(cbbPalette, length.out = length(unique(emb$celltype))),
                    unique(emb$celltype))) +
                theme_bw() + theme(aspect.ratio = 1) +
                labs(title = pipe, colour = "cell type")
        })
        if (length(panels) > 0) {
            p_umap = wrap_plots(panels, nrow = 1) + plot_layout(guides = "collect")
            pp_save_pdf(p_umap, pdir, "bio_umap_celltype",
                        width = 3.2 * length(panels), height = 3.2)
        }
    }

    ## Per-pipeline runtime and peak RSS from benchmark files.
    bench_dir = file.path(wd, "benchmarks")
    if (dir.exists(bench_dir)) {
        pipe_time = load_pipeline_time(bench_dir, aligners)
        pipe_mem = load_pipeline_memory(bench_dir, aligners)
        pp_save_csv(pipe_time, pdir, "bio_perf_time")
        pp_save_csv(pipe_mem, pdir, "bio_perf_memory")
        p_t = ggplot(pipe_time, aes(pipeline, total_min, fill = pipeline)) +
            geom_col() + scale_fill_manual(values = aligner_colours) +
            theme_bw() + theme(legend.position = "none") +
            labs(x = NULL, y = "wall-clock (min)")
        p_m = ggplot(pipe_mem, aes(pipeline, peak_rss_gb, fill = pipeline)) +
            geom_col() + scale_fill_manual(values = aligner_colours) +
            theme_bw() + theme(legend.position = "none") +
            labs(x = NULL, y = "peak RSS (GB)")
        pp_save_pdf(p_t + p_m, pdir, "bio_perf", width = 6, height = 3)
    }

    message("biology materials written to ", pdir)
}

## Benchmarks kind. Shared pipeline-agnostic plots and tables.
run_benchmarks = function(opt) {
    wd = opt$working_dir
    aligners = strsplit(opt$aligners, ",")[[1]]
    if (isTRUE(as.logical(opt$has_sbg))) aligners = c(aligners, "sbg")
    pdir = setup_paper_dir(wd, NULL)

    bench_dir = file.path(wd, "benchmarks")
    if (!dir.exists(bench_dir)) stop("benchmark dir missing: ", bench_dir)
    bm = load_benchmarks(bench_dir, aligners)
    pp_save_csv(bm, pdir, "bench_full")

    pipe_time = bm[pipeline %in% aligners,
                   .(total_min = sum(minutes)), by = pipeline]
    pipe_mem = bm[pipeline %in% aligners,
                  .(peak_rss_gb = max(max_rss_gb, na.rm = TRUE)), by = pipeline]
    pp_save_csv(pipe_time, pdir, "bench_total_time")
    pp_save_csv(pipe_mem, pdir, "bench_peak_memory")

    p_t = ggplot(pipe_time, aes(pipeline, total_min, fill = pipeline)) +
        geom_col() + geom_text(aes(label = round(total_min, 1)),
                               vjust = -0.3, size = 3) +
        scale_fill_manual(values = aligner_colours) +
        theme_bw() + theme(legend.position = "none") +
        labs(x = NULL, y = "total time (min)")
    p_m = ggplot(pipe_mem, aes(pipeline, peak_rss_gb, fill = pipeline)) +
        geom_col() + geom_text(aes(label = round(peak_rss_gb, 1)),
                               vjust = -0.3, size = 3) +
        scale_fill_manual(values = aligner_colours) +
        theme_bw() + theme(legend.position = "none") +
        labs(x = NULL, y = "peak RSS (GB)")
    pp_save_pdf(p_t + p_m, pdir, "bench_total", width = 6, height = 3)

    bm_steps = bm[pipeline %in% aligners]
    if (nrow(bm_steps) > 0) {
        p_step = ggplot(bm_steps, aes(reorder(step, minutes), minutes,
                                      fill = pipeline)) +
            geom_col(width = 0.7) + coord_flip() +
            scale_fill_manual(values = aligner_colours) +
            theme_bw() +
            labs(x = NULL, y = "wall-clock time (min)", fill = NULL)
        pp_save_pdf(p_step, pdir, "bench_step", width = 6, height = 4.5)
    }

    message("benchmark materials written to ", pdir)
}

## Benchmark helpers, shared across kinds.
load_benchmarks = function(bench_dir, aligners) {
    fs = list.files(bench_dir, pattern = "\\.txt$", full.names = TRUE)
    bm = rbindlist(Filter(Negate(is.null), lapply(fs, function(fn) {
        tryCatch({ dt = fread(fn, sep = "\t"); dt[, file := basename(fn)]; dt },
                 error = function(e) NULL)
    })), fill = TRUE)
    if (nrow(bm) == 0) return(bm)
    bm[, step := gsub("\\.txt$", "", file)]
    bm[, step := sub("_[A-Za-z0-9]+$", "", step)]
    classify = function(s) {
        if (grepl("starsolo|^star_", s)) return("starsolo")
        if (grepl("kallisto|bustools", s)) return("kallisto")
        if (grepl("alevin", s)) return("alevin")
        if (grepl("sbg", s)) return("sbg")
        "shared"
    }
    bm[, pipeline := vapply(step, classify, character(1))]
    bm[, minutes := s / 60]
    bm[, max_rss_gb := max_rss / 1024]
    bm[]
}

load_pipeline_time = function(bench_dir, aligners) {
    bm = load_benchmarks(bench_dir, aligners)
    bm[pipeline %in% aligners, .(total_min = sum(minutes)), by = pipeline]
}

load_pipeline_memory = function(bench_dir, aligners) {
    bm = load_benchmarks(bench_dir, aligners)
    bm[pipeline %in% aligners,
       .(peak_rss_gb = max(max_rss_gb, na.rm = TRUE)), by = pipeline]
}

main = function() {
    opt = parse_args()
    switch(opt$kind,
        simulation = run_simulation(opt),
        biology    = run_biology(opt),
        benchmarks = run_benchmarks(opt))
}

if (!interactive()) main()
