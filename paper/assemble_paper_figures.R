#!/usr/bin/env Rscript
## assemble_paper_figures.R
##
## Reads the RDS outputs produced by the main rhapsodist pipeline and writes
## manuscript-ready CSV and PDF artefacts under <working_dir>/paper/<sample>/
## (or <working_dir>/paper/ for benchmark-only outputs).
##
## Usage:
##   Rscript paper/assemble_paper_figures.R --kind simulation \
##       --working_dir output/simul --sample simulated \
##       --aligners starsolo,kallisto,alevin,sbg --n_expected_cells 1000
##
##   Rscript paper/assemble_paper_figures.R --kind biology \
##       --working_dir output/sendoel2024 --sample sample_16_wta_p60 \
##       --aligners starsolo,kallisto,alevin
##
##   Rscript paper/assemble_paper_figures.R --kind benchmarks \
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
                   choices = c("simulation", "biology", "benchmarks",
                               "linker_qc"))
    p$add_argument("--working_dir", required = TRUE)
    p$add_argument("--sample", default = "")
    p$add_argument("--aligners", default = "starsolo,kallisto,alevin")
    p$add_argument("--has_sbg", default = "false")
    p$add_argument("--n_expected_cells", type = "integer", default = 0L)
    p$add_argument("--markers_file", default = "")
    p$add_argument("--bench_prefix", default = "fig_sim")
    p$add_argument("--use_case",
                   choices = c("sendoel", "hela"),
                   default = "sendoel")
    p$parse_args()
}

aligner_colours = c(starsolo = "#E69F00", kallisto = "#56B4E9",
                    alevin = "#009E73", sbg = "#CC79A7",
                    truth = "#000000")

## Shared paper theme. Uses cowplot::theme_cowplot when available (no grids,
## harmonised fonts, paper-friendly), otherwise a theme_bw fallback with the
## background grid removed.
paper_theme = function(base_size = 10) {
    if (requireNamespace("cowplot", quietly = TRUE)) {
        cowplot::theme_cowplot(font_size = base_size) +
            theme(plot.tag = element_text(face = "bold", size = base_size + 2),
                  legend.key.size = grid::unit(0.35, "cm"),
                  plot.margin = grid::unit(c(8, 12, 8, 12), "pt"))
    } else {
        theme_bw(base_size = base_size) +
            theme(panel.grid = element_blank(),
                  plot.tag = element_text(face = "bold", size = base_size + 2),
                  legend.key.size = grid::unit(0.35, "cm"),
                  plot.margin = grid::unit(c(8, 12, 8, 12), "pt"))
    }
}

## Shared text sizes applied to every panel via patchwork's `&` operator so
## axis, legend, and strip labels are consistent across figures.
paper_shared_theme = theme(axis.title  = element_text(size = 9),
                           axis.text   = element_text(size = 8),
                           legend.title = element_text(size = 9),
                           legend.text  = element_text(size = 8),
                           strip.text   = element_text(size = 9),
                           plot.subtitle = element_text(size = 8))

## Render an UpSetR plot to a temporary PNG, then return a real ggplot that
## displays it via annotation_raster. Wrapping the upset as a grid grob via
## wrap_elements breaks patchwork composition (other panels collapse), so the
## raster trick is the workaround.
upset_panel = function(set_list, width_in = 6, height_in = 4, dpi = 200,
                       title = NULL) {
    if (!requireNamespace("UpSetR", quietly = TRUE)) return(NULL)
    if (!requireNamespace("png", quietly = TRUE)) return(NULL)
    if (length(set_list) < 2) return(NULL)
    tmp = tempfile(fileext = ".png")
    grDevices::png(tmp, width = width_in * dpi, height = height_in * dpi,
                   res = dpi, bg = "white")
    print(UpSetR::upset(UpSetR::fromList(set_list), order.by = "freq",
                        nsets = length(set_list), text.scale = 1.1,
                        mb.ratio = c(0.55, 0.45),
                        point.size = 2.2, line.size = 0.7))
    grDevices::dev.off()
    ## Trim the surrounding white margin that UpSetR leaves around its
    ## grid so the panel occupies its patchwork cell instead of looking
    ## tiny. Prefer magick::image_trim when available; otherwise do a
    ## pure-R bounding-box crop on the raw PNG so the trim still happens
    ## without the optional system dependency.
    if (requireNamespace("magick", quietly = TRUE)) {
        ## fuzz="5%" catches near-white pixels (antialiasing halo around
        ## UpSetR's lines and text) that plain image_trim would leave
        ## behind as margin.
        im = magick::image_trim(magick::image_read(tmp), fuzz = 5)
        tmp2 = tempfile(fileext = ".png")
        magick::image_write(im, path = tmp2, format = "png")
        img = png::readPNG(tmp2)
    } else {
        img = png::readPNG(tmp)
        ## More aggressive threshold (0.97 vs 0.995) to clip near-white
        ## halo in the pure-R fallback.
        is_white = if (length(dim(img)) == 3) {
            apply(img[, , seq_len(min(3, dim(img)[3])), drop = FALSE],
                  c(1, 2), min) >= 0.97
        } else {
            img >= 0.97
        }
        non_white_rows = which(!apply(is_white, 1, all))
        non_white_cols = which(!apply(is_white, 2, all))
        if (length(non_white_rows) >= 2 && length(non_white_cols) >= 2) {
            img = img[min(non_white_rows):max(non_white_rows),
                      min(non_white_cols):max(non_white_cols), ,
                      drop = FALSE]
        }
    }
    plot_h = nrow(img); plot_w = ncol(img)
    p = ggplot() +
        annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                          interpolate = TRUE) +
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1),
                        expand = FALSE, clip = "off") +
        theme_void() +
        theme(aspect.ratio = plot_h / plot_w,
              plot.title = element_text(size = 10, face = "plain",
                                        hjust = 0.5,
                                        margin = margin(b = 2)),
              plot.margin = grid::unit(c(0, 0, 0, 0), "pt"))
    if (!is.null(title)) p = p + ggtitle(title)
    p
}

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
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
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
                pdir, "sim_bc_upset", width = 6, height = 4,
                keep_last_page = TRUE)
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

        p_mard = ggplot(mard_dt, aes(pipeline1, pipeline2, fill = mard)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.3f", mard)), size = 3) +
            scale_fill_viridis_c(option = "viridis", direction = -1,
                                 limits = c(0, NA),
                                 alpha = 0.75) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
            labs(x = "pipeline", y = "pipeline", fill = "MARD",
                 title = "Pseudobulk MARD on raw counts")
        pp_save_pdf(p_mard, pdir, "sim_pseudobulk_mard", width = 5, height = 4)

        cor_mat = cor(pb_mat, method = "pearson")
        cor_dt = as.data.table(as.table(cor_mat))
        setnames(cor_dt, c("pipeline1", "pipeline2", "pearson_r"))
        pp_save_csv(cor_dt, pdir, "sim_pseudobulk_correlation")
        p_pcor = ggplot(cor_dt, aes(pipeline1, pipeline2, fill = pearson_r)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.3f", pearson_r)), size = 3) +
            scale_fill_viridis_c(option = "viridis", direction = 1,
                                 limits = c(min(cor_dt$pearson_r), 1),
                                 alpha = 0.75) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
            labs(x = "pipeline", y = "pipeline", fill = "Pearson r",
                 title = "Pseudobulk Pearson r on raw counts")
        pp_save_pdf(p_pcor, pdir, "sim_pseudobulk_correlation",
                    width = 5, height = 4)

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
        labs(x = "log10(UMI + 1)", y = "density", colour = "pipeline")
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

    ## Sample-tag demultiplexing: aggregate per-read assignments, then
    ## compare the majority-vote call per cell against the simulator's
    ## ground-truth assignment (simulate/sampletag_assignments.txt).
    st_fn = file.path(wd, "sampletags", samp, "sampletag_counts.tsv.gz")
    truth_st_fn = file.path(wd, "simulate", "sampletag_assignments.txt")
    if (file.exists(st_fn) && file.size(st_fn) > 0) {
        st = tryCatch(fread(cmd = paste("zcat", shQuote(st_fn)), header = FALSE,
                            col.names = c("barcode", "umi", "tag", "mismatches")),
                      error = function(e) NULL)
        if (is.null(st) || nrow(st) == 0) st = NULL
    } else {
        st = NULL
    }
    if (!is.null(st)) {
        tag_reads = st[, .(reads = .N,
                           cells = uniqueN(barcode)), by = tag]
        setorder(tag_reads, -reads)
        pp_save_csv(tag_reads, pdir, "sim_sampletag_counts")
        tag_long = melt(tag_reads, id.vars = "tag",
                        variable.name = "metric", value.name = "n")
        p_st = ggplot(tag_long, aes(tag, n, fill = tag)) +
            geom_col() + geom_text(aes(label = n), vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
            facet_wrap(~metric, scales = "free_y") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  legend.position = "none") +
            labs(x = NULL, y = NULL)
        pp_save_pdf(p_st, pdir, "sim_sampletag_counts",
                    width = 5.5, height = 3)

        if (file.exists(truth_st_fn)) {
            truth_st = fread(truth_st_fn)
            setnames(truth_st, c("barcode", "truth_tag"))
            calls = st[, .N, by = .(barcode, tag)]
            setorder(calls, barcode, -N)
            pred = calls[, .(pred_tag = tag[1], pred_reads = N[1],
                             total_reads = sum(N)), by = barcode]
            pred[, pred_purity := pred_reads / total_reads]
            eval = merge(truth_st, pred, by = "barcode", all.x = TRUE)
            eval[is.na(pred_tag), pred_tag := "unassigned"]
            conf = eval[, .N, by = .(truth_tag, pred_tag)]
            pp_save_csv(conf, pdir, "sim_sampletag_confusion")
            acc = eval[, .(accuracy = mean(truth_tag == pred_tag,
                                           na.rm = TRUE),
                           n_cells = .N,
                           mean_purity = mean(pred_purity, na.rm = TRUE))]
            pp_save_csv(acc, pdir, "sim_sampletag_accuracy")
            p_conf = ggplot(conf, aes(pred_tag, truth_tag, fill = N)) +
                geom_tile(colour = "white") +
                geom_text(aes(label = N), size = 3) +
                scale_fill_distiller(palette = "Blues", direction = 1,
                                     limits = c(0, max(conf$N))) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
                labs(x = "predicted tag", y = "true tag", fill = "cells",
                     subtitle = sprintf("accuracy %.1f%% (n=%d)",
                                        100 * acc$accuracy, acc$n_cells))
            pp_save_pdf(p_conf, pdir, "sim_sampletag_confusion",
                        width = 4.5, height = 3.5)
        }
    }

    panel_theme = paper_theme(10)
    panels = list()
    if (exists("bc_lists", inherits = FALSE)) {
        panels$B = upset_panel(bc_lists, width_in = 12, height_in = 6)
    }
    if (exists("p_mard", inherits = FALSE)) {
        panels$D = p_mard + panel_theme +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  aspect.ratio = 1)
    }
    if (exists("p_conf", inherits = FALSE)) {
        panels$C = p_conf + panel_theme +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  plot.subtitle = element_text(size = 8),
                  aspect.ratio = 1)
    }
    if (exists("p_pcor", inherits = FALSE)) {
        panels$G = p_pcor + panel_theme +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  aspect.ratio = 1)
    }
    bench_dir = file.path(wd, "benchmarks")
    if (dir.exists(bench_dir)) {
        bm = load_benchmarks(bench_dir, aligners)
        if (nrow(bm) > 0) {
            pipe_time = bm[pipeline %in% aligners & !is_install_rule(file),
                           .(total_min = sum(minutes)), by = pipeline]
            pipe_mem = bm[pipeline %in% aligners & !is_install_rule(file),
                          .(peak_rss_gb = max(max_rss_gb, na.rm = TRUE)),
                          by = pipeline]
            bar_theme = panel_theme +
                theme(legend.position = "none",
                      plot.margin = grid::unit(c(2, 6, 2, 2), "pt"),
                      axis.title.y = element_text(margin = margin(r = 2)),
                      axis.text.x = element_text(angle = 30, hjust = 1))
            panels$E = ggplot(pipe_time, aes(pipeline, total_min,
                                             fill = pipeline)) +
                geom_col(width = 0.55) +
                geom_text(aes(label = sprintf("%.2f", total_min)),
                          vjust = -0.3, size = 3) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
                scale_fill_manual(values = aligner_colours) +
                bar_theme + labs(x = "pipeline", y = "total time (min)")
            panels$F = ggplot(pipe_mem, aes(pipeline, peak_rss_gb,
                                            fill = pipeline)) +
                geom_col(width = 0.55) +
                geom_text(aes(label = sprintf("%.2f", peak_rss_gb)),
                          vjust = -0.3, size = 3) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
                scale_fill_manual(values = aligner_colours) +
                bar_theme + labs(x = "pipeline", y = "peak RSS (GB)")
        }
    }
    blank = function() ggplot() + theme_void()
    for (k in c("B", "C", "D", "E", "F", "G")) {
        if (is.null(panels[[k]])) {
            warning("fig1 panel ", k, " missing; rendering blank. wd=", wd)
            panels[[k]] = blank()
        }
    }
    ## fig 1: narrative-ordered B-F (A = workflow schematic is overlaid
    ## manually on the composed PDF). The simulated pseudobulk Pearson r
    ## heatmap (panels$G) is kept as a standalone supplementary PDF and
    ## not composed into fig 1, since MARD in panel D already conveys
    ## pseudobulk agreement on the simulated data.
    fig1_design = paste("BBBCCCDDD",
                        "BBBCCCDDD",
                        "BBBCCCDDD",
                        "EEEEEFFFF",
                        "EEEEEFFFF", sep = "\n")
    fig1 = patchwork::wrap_plots(B = panels$B, C = panels$C, D = panels$D,
                                 E = panels$E, F = panels$F,
                                 design = fig1_design) +
        patchwork::plot_annotation(tag_levels = list(c("B", "C", "D",
                                                       "E", "F"))) &
        paper_shared_theme
    pp_save_pdf(fig1, pdir, "fig1_simulations_panels",
                width = 10.5, height = 7)

    message("simulation materials written to ", pdir)
}

## Biology kind. Consumes derived biology_*.rds files written by the
## production report and Seurat caches for UMAP plots.
## Canonical plot order: sbg first (if present), then alevin, kallisto, starsolo.
## Used everywhere pipeline drives axis, fill, or facet order so figures read
## consistently across samples.
canonical_pipeline_order = c("sbg", "alevin", "kallisto", "starsolo", "truth")

order_pipelines = function(x, present = NULL) {
    if (is.null(present)) present = unique(as.character(x))
    lv = intersect(canonical_pipeline_order, present)
    factor(x, levels = lv)
}

## Disjoint hue ranges per aligner for cluster UMAPs. Hue windows are
## centred on each aligner's Okabe-Ito colour (aligner_colours above),
## so cluster hues in G sit in the same family as the aligner's bars
## and violins elsewhere and never overlap across aligners.
##   starsolo  #E69F00  → hue ~45   (orange/yellow)
##   kallisto  #56B4E9  → hue ~210  (sky blue)
##   alevin    #009E73  → hue ~160  (bluish green)
##   sbg       #CC79A7  → hue ~345  (reddish pink)
cluster_hue_range = list(
    starsolo = c(30,  70),
    kallisto = c(190, 240),
    alevin   = c(130, 180),
    sbg      = c(325, 365))

cluster_palette_for = function(pipe, n_levels) {
    h_range = cluster_hue_range[[pipe]]
    if (is.null(h_range)) h_range = c(0, 360)
    if (n_levels <= 1) {
        return(hcl(h = mean(h_range), c = 90, l = 55))
    }
    ## Sweep hue, luminance, and chroma together so clusters within a
    ## single aligner's zone differ in tint/shade as well as hue.
    hues = seq(h_range[1], h_range[2], length.out = n_levels)
    lums = seq(40, 75, length.out = n_levels)
    chromas = seq(100, 70, length.out = n_levels)
    hcl(h = hues, c = chromas, l = lums)
}

run_biology = function(opt) {
    wd = opt$working_dir
    samp = opt$sample
    aligners = strsplit(opt$aligners, ",")[[1]]
    if (isTRUE(as.logical(opt$has_sbg))) aligners = c(aligners, "sbg")
    aligners = intersect(canonical_pipeline_order, aligners)
    pdir = setup_paper_dir(wd, samp)
    use_case = if (is.null(opt$use_case)) "sendoel" else opt$use_case
    stopifnot(use_case %in% c("sendoel", "hela"))
    is_hela = identical(use_case, "hela")
    fig_stem = switch(use_case, hela = "fig2_hela", sendoel = "fig3_sendoel")

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
            geom_col() +
            geom_text(aes(label = signif(value, 3)),
                      vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
            facet_wrap(~metric, scales = "free_y") +
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
        mard_dt[, pipeline1 := order_pipelines(pipeline1, aligners)]
        mard_dt[, pipeline2 := order_pipelines(pipeline2, aligners)]
        p = ggplot(mard_dt, aes(pipeline1, pipeline2, fill = mard * 100)) +
            geom_tile(colour = "white") +
            geom_text(aes(label = sprintf("%.1f", mard * 100)), size = 2.8) +
            scale_fill_viridis_c(option = "viridis", direction = -1,
                                 limits = c(0, 100), alpha = 0.75) +
            scale_y_discrete(limits = rev) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  aspect.ratio = 1) +
            labs(x = NULL, y = NULL, fill = "MARD (%)",
                 title = "Pseudobulk MARD")
        pp_save_pdf(p, pdir, "bio_pseudobulk_mard", width = 4.5, height = 3.5)

        ## Pairwise pseudobulk scatter, log10 counts, marker genes highlighted.
        markers = character(0)
        if (nzchar(opt$markers_file) && file.exists(opt$markers_file)) {
            marker_meta = read.table(opt$markers_file, header = TRUE, sep = "\t",
                                     stringsAsFactors = FALSE)
            markers = marker_meta$marker
        }
        if (ncol(pb_mat) >= 2) {
            pb_log = log10(pb_mat + 1)
            pb_dt = as.data.table(pb_log, keep.rownames = "gene")
            pb_dt[, is_marker := gene %in% markers]
            cor_pairs = combn(colnames(pb_mat), 2, simplify = FALSE)
            cor_rows = list()
            scatter_plots = lapply(cor_pairs, function(pair) {
                d = pb_dt[, .(gene, x = get(pair[1]), y = get(pair[2]), is_marker)]
                r = cor(d$x, d$y)
                cor_rows[[paste(pair, collapse = "_vs_")]] <<- data.table(
                    pipeline1 = pair[1], pipeline2 = pair[2],
                    n_genes = nrow(d), pearson_r = r)
                ggplot(d, aes(x = x, y = y)) +
                    pp_rasterise(geom_point(data = d[is_marker == FALSE],
                                            alpha = 0.15, size = 0.3,
                                            colour = "grey50")) +
                    geom_point(data = d[is_marker == TRUE],
                               colour = "#E69F00", size = 1.2) +
                    ggrepel::geom_text_repel(data = d[is_marker == TRUE],
                                             aes(label = gene), size = 2.4,
                                             colour = "#E69F00",
                                             max.overlaps = Inf,
                                             min.segment.length = 0,
                                             segment.size = 0.2,
                                             box.padding = 0.3) +
                    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                                colour = "grey30") +
                    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
                             label = sprintf("R = %.4f", r), size = 3) +
                    theme_bw() + theme(aspect.ratio = 1) +
                    labs(x = paste0(pair[1], " log10(count + 1)"),
                         y = paste0(pair[2], " log10(count + 1)"),
                         title = paste(pair[2], "vs", pair[1]))
            })
            pp_save_csv(rbindlist(cor_rows), pdir, "bio_pseudobulk_correlation")
            pp_save_pdf(wrap_plots(scatter_plots, nrow = 1), pdir,
                        "bio_pseudobulk_correlation",
                        width = 3.2 * length(scatter_plots), height = 3.2)

            ## Symmetric Pearson r heatmap, mirroring the simulation panel.
            bcor_mat = cor(pb_mat, method = "pearson")
            bcor_dt = as.data.table(as.table(bcor_mat))
            setnames(bcor_dt, c("pipeline1", "pipeline2", "pearson_r"))
            pp_save_csv(bcor_dt, pdir, "bio_pseudobulk_correlation_matrix")
            p_bcor = ggplot(bcor_dt,
                            aes(pipeline1, pipeline2, fill = pearson_r)) +
                geom_tile(colour = "white") +
                geom_text(aes(label = sprintf("%.3f", pearson_r)), size = 3) +
                scale_fill_viridis_c(option = "viridis", direction = 1,
                                     limits = c(0, 1),
                                     alpha = 0.75) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
                labs(x = "pipeline", y = "pipeline", fill = "Pearson r",
                     title = "Pseudobulk Pearson r on log-normalized counts")
            pp_save_pdf(p_bcor, pdir, "bio_pseudobulk_correlation_heatmap",
                        width = 5, height = 4)
        }
    }

    ## Barcode overlap across pipelines.
    if (file.exists(biords("cb_umi"))) {
        bc_dt = as.data.table(readRDS(biords("cb_umi")))
        bc_lists = split(bc_dt$barcode, bc_dt$pipeline)
        if (requireNamespace("UpSetR", quietly = TRUE) && length(bc_lists) >= 2) {
            pp_save_base_pdf(
                quote(print(UpSetR::upset(UpSetR::fromList(bc_lists),
                                          order.by = "freq",
                                          nsets = length(bc_lists),
                                          text.scale = 1.2))),
                pdir, "bio_bc_upset", width = 6, height = 4,
                keep_last_page = TRUE)
            set_sizes = data.table(set = names(bc_lists),
                                   size = vapply(bc_lists, length, integer(1)))
            pp_save_csv(set_sizes, pdir, "bio_bc_upset_sizes")
        }
    }

    ## Cell recovery bar plot (separate panel for parity with the simulation).
    if (file.exists(biords("qc"))) {
        qc_df = as.data.table(readRDS(biords("qc")))
        if ("n_cells" %in% colnames(qc_df)) {
            p_rec = ggplot(qc_df, aes(pipeline, n_cells, fill = pipeline)) +
                geom_col() +
                geom_text(aes(label = n_cells), vjust = -0.3, size = 3) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
                scale_fill_manual(values = aligner_colours) +
                theme_bw() + theme(legend.position = "none") +
                labs(x = NULL, y = "cells recovered")
            pp_save_pdf(p_rec, pdir, "bio_cell_counts", width = 4, height = 3)
            pp_save_csv(qc_df[, .(pipeline, n_cells)], pdir, "bio_cell_counts")
        }
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
                    labs(x = sprintf("mean log10(UMI+1): (%s + %s) / 2",
                                     pair[1], pair[2]),
                         y = sprintf("log10(UMI+1): %s - %s",
                                     pair[1], pair[2]),
                         title = paste(pair[1], "-", pair[2]),
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
                geom_text(aes(y = ci_hi, label = sprintf("%.3f", ari)),
                          vjust = -0.4, size = 3) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
                labs(x = NULL, y = ylab)
            pp_save_pdf(p, pdir, name, width = 4.5, height = 3.2)
        }
        ari_bar(cluster_ari, "bio_cluster_ari", "cluster ARI")
        ari_bar(ct_ari,      "bio_celltype_ari", "cell-type ARI")

        ## Combined cluster + cell-type ARI bar panel. Used as fig 3 panel F
        ## on the sendoel use case so the main figure exposes ARI directly
        ## next to the UMAPs and performance bars. Only built when both
        ## tables are present (i.e. the marker-voting celltype exists, which
        ## is the case for sendoel but not hela).
        p_ari_combo = NULL
        if (!is.null(cluster_ari) && !is.null(ct_ari) &&
            nrow(cluster_ari) > 0 && nrow(ct_ari) > 0) {
            ari_combo = rbind(
                cluster_ari[, .(pair = paste(pipeline1, pipeline2,
                                             sep = " vs "),
                                ari, ci_lo, ci_hi, kind = "cluster")],
                ct_ari[,      .(pair = paste(pipeline1, pipeline2,
                                             sep = " vs "),
                                ari, ci_lo, ci_hi, kind = "cell type")])
            pp_save_csv(ari_combo, pdir, "bio_ari_combo")
            p_ari_combo = ggplot(ari_combo,
                                 aes(pair, ari, fill = kind)) +
                geom_col(position = position_dodge(0.75), width = 0.65) +
                geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                              position = position_dodge(0.75),
                              width = 0.25, linewidth = 0.4) +
                geom_text(aes(y = ci_hi, label = sprintf("%.3f", ari)),
                          position = position_dodge(0.75),
                          vjust = -0.4, size = 2.8) +
                scale_y_continuous(expand = expansion(mult = c(0.05, 0.12)),
                                   limits = c(0, 1)) +
                scale_fill_manual(values = c(cluster = "#0072B2",
                                             `cell type` = "#E69F00")) +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 30, hjust = 1),
                      legend.position = "top") +
                labs(x = NULL, y = "ARI", fill = NULL)
            pp_save_pdf(p_ari_combo, pdir, "bio_ari_combo",
                        width = 5.5, height = 3.2)
        }

        ## Pairwise confusion-matrix heatmaps underlying each ARI. Rows are
        ## labels from pipeline A, columns from pipeline B, cell entries are
        ## cell counts on shared barcodes.
        draw_confusion = function(key, name, xlab) {
            pipes = sort(unique(clusters_dt$pipeline))
            if (length(pipes) < 2) return(invisible(NULL))
            pairs = combn(pipes, 2, simplify = FALSE)
            panels = lapply(pairs, function(pair) {
                d1 = clusters_dt[pipeline == pair[1], .(barcode, a = get(key))]
                d2 = clusters_dt[pipeline == pair[2], .(barcode, b = get(key))]
                sh = merge(d1, d2, by = "barcode")
                if (nrow(sh) < 10) return(NULL)
                tab = as.data.table(sh[, .N, by = .(a, b)])
                ggplot(tab, aes(a, b, fill = N)) +
                    geom_tile(colour = "white") +
                    scale_fill_distiller(palette = "Blues", direction = 1) +
                    theme_bw(base_size = 11) +
                    theme(aspect.ratio = 1,
                          axis.text.x = element_text(angle = 45,
                                                     hjust = 1, vjust = 1),
                          legend.position = "right",
                          plot.margin = margin(4, 4, 4, 4)) +
                    labs(x = pair[1], y = pair[2], fill = "cells",
                         title = paste(pair[2], "vs", pair[1]))
            })
            panels = Filter(Negate(is.null), panels)
            if (length(panels) == 0) return(invisible(NULL))
            p = wrap_plots(panels, nrow = 1)
            pp_save_pdf(p, pdir, name,
                        width = 14, height = 14 / length(panels))
        }
        draw_confusion("cluster_prefixed", "bio_cluster_confusion", "cluster")
        draw_confusion("celltype",         "bio_celltype_confusion", "cell type")
    }

    ## UMAP panels. For the sendoel use case, D is celltype UMAP and E is
    ## cluster UMAP. For the hela use case (no marker celltype on a clonal
    ## line) D is an ARI heatmap of cluster concordance across aligners and
    ## E is a UMAP coloured by Seurat cell cycle phase. The Seurat caches
    ## are saved by 04_biology.Rmd; the per-barcode table with cluster,
    ## celltype, and cell cycle phase columns is in the clusters RDS.
    seu_fns = setNames(file.path(wd, sprintf("%s_biology_%s_seurat.rds", samp,
                                             aligners)), aligners)
    clusters_fn = biords("clusters")
    if (all(file.exists(seu_fns)) && file.exists(clusters_fn) &&
        requireNamespace("Seurat", quietly = TRUE)) {
        cbbPalette = c("#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        clusters_dt = as.data.table(readRDS(clusters_fn))

        if (!is_hela) {
            ct_levels = c(sort(setdiff(unique(clusters_dt$celltype), "none")),
                          "none")
            ct_colours = setNames(
                c(rep(cbbPalette, length.out = length(ct_levels) - 1), "grey85"),
                ct_levels)
            panels = lapply(names(seu_fns), function(pipe) {
                so = readRDS(seu_fns[[pipe]])
                emb = as.data.frame(Seurat::Embeddings(so, "umap"))
                colnames(emb) = c("UMAP_1", "UMAP_2")
                emb$barcode = rownames(emb)
                ct_map = clusters_dt[pipeline == pipe,
                                     setNames(celltype, barcode)]
                emb$celltype = factor(
                    ifelse(emb$barcode %in% names(ct_map),
                           ct_map[emb$barcode], "none"),
                    levels = ct_levels)
                ggplot(emb, aes(UMAP_1, UMAP_2, colour = celltype)) +
                    pp_rasterise(geom_point(size = 0.3, alpha = 0.8)) +
                    scale_colour_manual(values = ct_colours, drop = FALSE) +
                    guides(colour = guide_legend(
                        override.aes = list(size = 2.5, alpha = 1))) +
                    theme_bw() + theme(aspect.ratio = 1,
                                       panel.grid = element_blank()) +
                    labs(title = pipe, colour = "cell type")
            })
            if (length(panels) > 0) {
                p_umap = wrap_plots(panels, nrow = 1) +
                    plot_layout(guides = "collect") +
                    plot_annotation(title = "    Cell embeddings by annotation",
                                    theme = theme(plot.title = element_text(
                                        size = 11, face = "plain",
                                        margin = margin(l = 30, b = 4))))
                pp_save_pdf(p_umap, pdir, "bio_umap_celltype",
                            width = 3.2 * length(panels), height = 3.2)
            }

            cl_panels = lapply(names(seu_fns), function(pipe) {
                so = readRDS(seu_fns[[pipe]])
                emb = as.data.frame(Seurat::Embeddings(so, "umap"))
                colnames(emb) = c("UMAP_1", "UMAP_2")
                emb$barcode = rownames(emb)
                cl_map = clusters_dt[pipeline == pipe,
                                     setNames(as.character(cluster_prefixed),
                                              barcode)]
                emb$cluster = factor(
                    ifelse(emb$barcode %in% names(cl_map),
                           cl_map[emb$barcode], NA))
                ggplot(emb, aes(UMAP_1, UMAP_2, colour = cluster)) +
                    pp_rasterise(geom_point(size = 0.3, alpha = 0.8)) +
                    guides(colour = guide_legend(ncol = 2,
                                                 override.aes = list(size = 1.5))) +
                    theme_bw() + theme(aspect.ratio = 1,
                                       panel.grid = element_blank(),
                                       legend.position = "right",
                                       legend.key.size = grid::unit(0.3, "cm"),
                                       legend.text = element_text(size = 7)) +
                    labs(title = pipe, colour = "Louvain")
            })
            if (length(cl_panels) > 0) {
                p_cl = wrap_plots(cl_panels, nrow = 1) +
                    plot_annotation(title = "    Cell embeddings by cluster",
                                    theme = theme(plot.title = element_text(
                                        size = 11, face = "plain",
                                        margin = margin(l = 30, t = 2, b = 6))))
                pp_save_pdf(p_cl, pdir, "bio_umap_cluster",
                            width = 3.6 * length(cl_panels), height = 3.2)
            }
        } else {
            ## hela use case. Fig 2 bottom block will stack two 1x4 UMAP rows
            ## (cluster on top, cell cycle phase below) on the left 8/12 cols
            ## and stack time/memory bars on the right 4/12 cols. The cluster
            ## and phase ARI heatmaps fill panels E and F in the middle row.

            ## Cluster ARI heatmap (panel E).
            if (requireNamespace("mclust", quietly = TRUE) &&
                length(aligners) >= 2) {
                pair_grid = expand.grid(a = aligners, b = aligners,
                                        stringsAsFactors = FALSE)
                ari_rows = rbindlist(lapply(seq_len(nrow(pair_grid)),
                                            function(i) {
                    aa = pair_grid$a[i]; bb = pair_grid$b[i]
                    d1 = clusters_dt[pipeline == aa,
                                     .(barcode, cl = cluster_prefixed)]
                    d2 = clusters_dt[pipeline == bb,
                                     .(barcode, cl = cluster_prefixed)]
                    sh = merge(d1, d2, by = "barcode",
                               suffixes = c(".a", ".b"))
                    ari = if (nrow(sh) < 10) NA_real_
                          else mclust::adjustedRandIndex(sh$cl.a, sh$cl.b)
                    data.table(pipeline1 = aa, pipeline2 = bb, ari = ari)
                }))
                pp_save_csv(ari_rows, pdir, "bio_cluster_ari_matrix")
                ari_rows[, pipeline1 := order_pipelines(pipeline1, aligners)]
                ari_rows[, pipeline2 := order_pipelines(pipeline2, aligners)]
                p_ari_heat = ggplot(ari_rows,
                                    aes(pipeline1, pipeline2, fill = ari)) +
                    geom_tile(colour = "white") +
                    geom_text(aes(label = ifelse(is.na(ari), "",
                                                 sprintf("%.3f", ari))),
                              size = 2.8) +
                    scale_fill_viridis_c(option = "viridis", direction = 1,
                                         limits = c(0, 1),
                                         alpha = 0.75, na.value = "grey90") +
                    scale_y_discrete(limits = rev) +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 30, hjust = 1),
                          aspect.ratio = 1) +
                    labs(x = NULL, y = NULL, fill = "ARI",
                         title = "Cluster ARI")
                pp_save_pdf(p_ari_heat, pdir, "bio_cluster_ari_matrix",
                            width = 5, height = 4)
            }

            ## Cell-cycle phase ARI heatmap (panel F).
            if (requireNamespace("mclust", quietly = TRUE) &&
                "Phase" %in% colnames(clusters_dt) &&
                any(!is.na(clusters_dt$Phase)) &&
                length(aligners) >= 2) {
                ph_grid = expand.grid(a = aligners, b = aligners,
                                      stringsAsFactors = FALSE)
                ph_rows = rbindlist(lapply(seq_len(nrow(ph_grid)),
                                           function(i) {
                    aa = ph_grid$a[i]; bb = ph_grid$b[i]
                    d1 = clusters_dt[pipeline == aa,
                                     .(barcode, ph = Phase)]
                    d2 = clusters_dt[pipeline == bb,
                                     .(barcode, ph = Phase)]
                    sh = merge(d1, d2, by = "barcode",
                               suffixes = c(".a", ".b"))
                    sh = sh[!is.na(ph.a) & !is.na(ph.b)]
                    ari = if (nrow(sh) < 10) NA_real_
                          else mclust::adjustedRandIndex(sh$ph.a, sh$ph.b)
                    data.table(pipeline1 = aa, pipeline2 = bb, ari = ari)
                }))
                pp_save_csv(ph_rows, pdir, "bio_phase_ari_matrix")
                ph_rows[, pipeline1 := order_pipelines(pipeline1, aligners)]
                ph_rows[, pipeline2 := order_pipelines(pipeline2, aligners)]
                p_phase_ari = ggplot(ph_rows,
                                     aes(pipeline1, pipeline2, fill = ari)) +
                    geom_tile(colour = "white") +
                    geom_text(aes(label = ifelse(is.na(ari), "",
                                                 sprintf("%.3f", ari))),
                              size = 2.8) +
                    scale_fill_viridis_c(option = "viridis", direction = 1,
                                         limits = c(0, 1),
                                         alpha = 0.75,
                                         na.value = "grey90") +
                    scale_y_discrete(limits = rev) +
                    theme_bw() +
                    theme(axis.text.x = element_text(angle = 30, hjust = 1),
                          aspect.ratio = 1) +
                    labs(x = NULL, y = NULL, fill = "ARI",
                         title = "Cellcycle phase ARI")
                pp_save_pdf(p_phase_ari, pdir, "bio_phase_ari_matrix",
                            width = 5, height = 4)
            }

            ## Cluster UMAPs (panel G): one per aligner, horizontal row.
            ## Per-aligner Louvain labels are distinct (s_0, k_0, ...), so
            ## legends cannot be collected; render one legend below each
            ## UMAP in two columns to match the UMAP width. Colors per
            ## aligner come from cluster_palette_for() which uses disjoint
            ## HCL hue zones so no cluster color repeats across aligners.
            ## Centroid labels help distinguish clusters within an aligner's
            ## narrow hue band. UMAP_1/UMAP_2 axis titles appear only on the
            ## leftmost panel.
            cl_panels_hela = lapply(seq_along(aligners), function(i) {
                pipe = aligners[i]
                if (!file.exists(seu_fns[[pipe]])) return(NULL)
                so = readRDS(seu_fns[[pipe]])
                emb = as.data.frame(Seurat::Embeddings(so, "umap"))
                colnames(emb) = c("UMAP_1", "UMAP_2")
                emb$barcode = rownames(emb)
                cl_map = clusters_dt[pipeline == pipe,
                                     setNames(as.character(cluster_prefixed),
                                              barcode)]
                emb$cluster = factor(
                    ifelse(emb$barcode %in% names(cl_map),
                           cl_map[emb$barcode], NA))
                pal = cluster_palette_for(pipe, nlevels(emb$cluster))
                emb_dt = as.data.table(emb)
                centroids = emb_dt[!is.na(cluster),
                                   .(UMAP_1 = median(UMAP_1),
                                     UMAP_2 = median(UMAP_2)),
                                   by = cluster]
                axis_title_theme = if (i == 1) element_text(size = 9)
                                   else element_blank()
                ggplot(emb, aes(UMAP_1, UMAP_2, colour = cluster)) +
                    pp_rasterise(geom_point(size = 0.3, alpha = 0.8)) +
                    scale_colour_manual(values = pal, na.value = "grey85") +
                    ggrepel::geom_text_repel(
                        data = centroids,
                        aes(label = cluster),
                        colour = "black",
                        size = 2.5,
                        fontface = "bold",
                        bg.colour = "white",
                        bg.r = 0.15,
                        box.padding = 0.25,
                        min.segment.length = 0,
                        segment.size = 0.2,
                        max.overlaps = Inf,
                        seed = 1L,
                        show.legend = FALSE) +
                    guides(colour = guide_legend(ncol = 2,
                                                 override.aes = list(size = 1.5))) +
                    theme_bw() +
                    theme(aspect.ratio = 1,
                          panel.grid = element_blank(),
                          axis.title.x = axis_title_theme,
                          axis.title.y = axis_title_theme,
                          legend.position = "bottom",
                          legend.key.size = grid::unit(0.3, "cm"),
                          legend.text = element_text(size = 7),
                          legend.title = element_text(size = 8)) +
                    labs(title = pipe, colour = "Louvain")
            })
            cl_panels_hela = Filter(Negate(is.null), cl_panels_hela)
            if (length(cl_panels_hela) > 0) {
                p_cluster_row = wrap_plots(cl_panels_hela, nrow = 1)
                pp_save_pdf(p_cluster_row, pdir, "bio_umap_cluster",
                            width = 3.6 * length(cl_panels_hela), height = 3.6)
            }

            ## Cell-cycle phase UMAPs (panel H): one per aligner, horizontal row.
            ## Shared categorical (G1, S, G2M) → collected legend, single row below.
            if ("Phase" %in% colnames(clusters_dt) &&
                any(!is.na(clusters_dt$Phase))) {
                phase_levels = c("G1", "S", "G2M")
                phase_colours = setNames(c("#1B9E77", "#D95F02", "#7570B3"),
                                         phase_levels)
                cc_panels = lapply(seq_along(aligners), function(i) {
                    pipe = aligners[i]
                    if (!file.exists(seu_fns[[pipe]])) return(NULL)
                    so = readRDS(seu_fns[[pipe]])
                    emb = as.data.frame(Seurat::Embeddings(so, "umap"))
                    colnames(emb) = c("UMAP_1", "UMAP_2")
                    emb$barcode = rownames(emb)
                    ph_map = clusters_dt[pipeline == pipe,
                                         setNames(as.character(Phase),
                                                  barcode)]
                    emb$phase = factor(
                        ifelse(emb$barcode %in% names(ph_map),
                               ph_map[emb$barcode], NA),
                        levels = phase_levels)
                    axis_title_theme = if (i == 1) element_text(size = 9)
                                       else element_blank()
                    ggplot(emb, aes(UMAP_1, UMAP_2, colour = phase)) +
                        pp_rasterise(geom_point(size = 0.3, alpha = 0.8)) +
                        scale_colour_manual(values = phase_colours,
                                            na.value = "grey85",
                                            drop = FALSE) +
                        guides(colour = guide_legend(
                            nrow = 1,
                            override.aes = list(size = 2.5, alpha = 1))) +
                        theme_bw() +
                        theme(aspect.ratio = 1,
                              panel.grid = element_blank(),
                              axis.title.x = axis_title_theme,
                              axis.title.y = axis_title_theme,
                              legend.key.size = grid::unit(0.3, "cm"),
                              legend.text = element_text(size = 7),
                              legend.title = element_text(size = 8)) +
                        labs(title = pipe, colour = "Cell cycle phase")
                })
                cc_panels = Filter(Negate(is.null), cc_panels)
                if (length(cc_panels) > 0) {
                    p_phase_row = wrap_plots(cc_panels, nrow = 1) +
                        plot_layout(guides = "collect") &
                        theme(legend.position = "bottom")
                    pp_save_pdf(p_phase_row, pdir, "bio_umap_phase",
                                width = 3.6 * length(cc_panels), height = 3.6)
                }
            }
        }
    }

    per_cell_qc = NULL
    p_per_cell_cor = NULL
    sce_fns = setNames(file.path(wd, aligners, samp,
                                 paste0(samp, "_", aligners, "_sce.rds")),
                       aligners)
    sce_list = NULL
    if (all(file.exists(sce_fns))) {
        sce_list = lapply(sce_fns, readRDS)
        per_cell_qc = rbindlist(lapply(names(sce_list), function(pipe) {
            sce = sce_list[[pipe]]
            counts_m = counts(sce)
            total_umi = colSums(counts_m)
            n_genes = colSums(counts_m > 0)
            gene_names = rownames(counts_m)
            mito_idx = grepl("^(mt-|MT-|Mt-)", gene_names)
            if (!any(mito_idx)) {
                rd = SummarizedExperiment::rowData(sce)
                sym_col = intersect(c("Symbol", "symbol", "gene_name",
                                      "gene_symbol", "name"), colnames(rd))
                if (length(sym_col) > 0) {
                    syms = as.character(rd[[sym_col[1]]])
                    mito_idx = grepl("^(mt-|MT-|Mt-)", syms)
                }
            }
            if (!any(mito_idx)) {
                rd = SummarizedExperiment::rowData(sce)
                chr_col = intersect(c("chr", "seqnames", "chromosome"),
                                    colnames(rd))
                if (length(chr_col) > 0) {
                    chr_vals = as.character(rd[[chr_col[1]]])
                    mito_idx = chr_vals %in% c("MT", "chrM", "M", "Mt", "chrMT")
                }
            }
            mito_pct = if (any(mito_idx))
                100 * colSums(counts_m[mito_idx, , drop = FALSE]) / pmax(total_umi, 1)
            else rep(NA_real_, ncol(sce))
            data.table(pipeline = pipe, total_umi = total_umi,
                       n_genes = n_genes, mito_pct = mito_pct)
        }))
        pp_save_csv(per_cell_qc, pdir, "bio_per_cell_qc")
    }

    ## Per-cell cross-aligner Pearson r on log1p counts over shared barcodes
    ## and shared genes. Mirrors the density plot in 02_comparison.Rmd, stored
    ## here as both a CSV summary and a standalone PDF so it can be composed
    ## into the biology figure panel. Shared-cell set is downsampled to at
    ## most max_cells_per_cell_cor per pair with a fixed seed so the density
    ## stays cheap to compute on large real datasets and is reproducible.
    ## Only computed for the HeLa use case because fig3 (sendoel) does not
    ## include this panel and the dense matrix materialisation is expensive
    ## on full experimental SCEs.
    max_cells_per_cell_cor = 500L
    if (is_hela && !is.null(sce_list) && length(sce_list) >= 2) {
        col_pearson = function(A, B) {
            Am = colMeans(A); Bm = colMeans(B)
            Ac = sweep(A, 2, Am, "-"); Bc = sweep(B, 2, Bm, "-")
            num = colSums(Ac * Bc)
            denom = sqrt(colSums(Ac^2) * colSums(Bc^2))
            ifelse(denom > 0, num / denom, NA_real_)
        }
        ## Harmonize gene identifiers across aligners so SBG (rownames are
        ## gene symbols) can be matched against starsolo/kallisto/alevin
        ## (rownames are Ensembl IDs; gene symbols live in rowData). We
        ## match on keys and subset the counts matrix by integer index,
        ## rather than overwriting rownames on the SCE wrapper, because the
        ## HDF5-backed assay does not always inherit the updated dimnames.
        ## Picking the rowData column: prefer explicitly-named symbol
        ## columns, otherwise the highest-cardinality column (skips static
        ## columns like "type"="Gene" or "value"="Expression" added by
        ## STARsolo's features.tsv import).
        gene_key = function(sce) {
            rd = SummarizedExperiment::rowData(sce)
            for (col in c("Symbol", "symbol", "gene_name", "gene_symbol")) {
                if (col %in% colnames(rd)) {
                    v = as.character(rd[[col]])
                    if (length(v) && all(nzchar(v))) return(v)
                }
            }
            n = nrow(sce)
            best_col = NA_character_
            best_u = 1L
            for (col in colnames(rd)) {
                v = as.character(rd[[col]])
                if (length(v) != n || !all(nzchar(v))) next
                u = length(unique(v))
                if (u > best_u && u >= max(10, n * 0.5)) {
                    best_col = col
                    best_u = u
                }
            }
            if (!is.na(best_col)) return(as.character(rd[[best_col]]))
            rn = rownames(sce)
            if (!is.null(rn) && length(rn) && all(nzchar(rn))) return(rn)
            as.character(seq_len(n))
        }
        gene_keys = lapply(sce_list, gene_key)
        pair_combos = combn(names(sce_list), 2, simplify = FALSE)
        per_cell_cor = rbindlist(lapply(pair_combos, function(pair) {
            a = pair[1]; b = pair[2]
            shared_cells = intersect(colnames(sce_list[[a]]),
                                     colnames(sce_list[[b]]))
            key_a = gene_keys[[a]]
            key_b = gene_keys[[b]]
            shared_genes = intersect(key_a, key_b)
            if (length(shared_cells) < 2 || length(shared_genes) < 2) return(NULL)
            ## integer row indices into each SCE's counts matrix; first-match
            ## wins on duplicates, which is fine for the density plot.
            ix_a = match(shared_genes, key_a)
            ix_b = match(shared_genes, key_b)
            n_total = length(shared_cells)
            if (n_total > max_cells_per_cell_cor) {
                shared_cells = shared_cells[order(shared_cells)]
                set.seed(1L)
                shared_cells = sort(sample(shared_cells,
                                           max_cells_per_cell_cor,
                                           replace = FALSE))
            }
            A = log1p(as.matrix(counts(sce_list[[a]])[ix_a,
                                                      shared_cells, drop = FALSE]))
            B = log1p(as.matrix(counts(sce_list[[b]])[ix_b,
                                                      shared_cells, drop = FALSE]))
            data.table(pair = paste(a, "vs", b),
                       cell = shared_cells,
                       correlation = col_pearson(A, B),
                       n_sampled = length(shared_cells),
                       n_shared_total = n_total)
        }), fill = TRUE)
        if (nrow(per_cell_cor) > 0) {
            valid_pc = per_cell_cor[!is.na(correlation)]
            pc_summary = valid_pc[, .(n_sampled = .N,
                                      n_shared_total = n_shared_total[1],
                                      median_r = median(correlation),
                                      mean_r   = mean(correlation),
                                      q25      = quantile(correlation, 0.25),
                                      q75      = quantile(correlation, 0.75)),
                                  by = pair]
            pp_save_csv(pc_summary, pdir, "bio_per_cell_correlation_summary")
            saveRDS(valid_pc, file.path(pdir, "bio_per_cell_correlation.rds"))
            ## Palette chosen to stay off the aligner hues (Okabe-Ito orange,
            ## sky-blue, bluish-green, pink). Up to 6 curves, one per aligner
            ## pair.
            pair_colours = c("#8B0000", "#4B0082", "#556B2F",
                             "#8B4513", "#2F4F4F", "#6A5ACD")
            p_per_cell_cor = ggplot(valid_pc,
                                    aes(x = correlation, colour = pair,
                                        fill = pair)) +
                geom_density(alpha = 0.25, linewidth = 0.7) +
                scale_colour_manual(values = pair_colours) +
                scale_fill_manual(values = pair_colours) +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      aspect.ratio = 1) +
                labs(x = "per-cell Pearson r", y = "density",
                     colour = NULL, fill = NULL,
                     title = "Per-cell cross-aligner correlation",
                     subtitle = sprintf(
                         "Pearson on log1p counts, shared genes, up to %d cells per pair",
                         max_cells_per_cell_cor))
            pp_save_pdf(p_per_cell_cor, pdir, "bio_per_cell_correlation",
                        width = 4.5, height = 4.5)
        }
    }
    rm(sce_list); invisible(gc(verbose = FALSE))

    ## Per-pipeline runtime and peak RSS from benchmark files.
    pipe_time = NULL; pipe_mem = NULL
    bench_dir = file.path(wd, "benchmarks")
    if (dir.exists(bench_dir)) {
        pipe_time = load_pipeline_time(bench_dir, aligners)
        pipe_mem = load_pipeline_memory(bench_dir, aligners)
        pp_save_csv(pipe_time, pdir, "bio_perf_time")
        pp_save_csv(pipe_mem, pdir, "bio_perf_memory")
        pipe_time[, pipeline := order_pipelines(pipeline, aligners)]
        pipe_mem[, pipeline := order_pipelines(pipeline, aligners)]
        p_t = ggplot(pipe_time, aes(pipeline, total_min, fill = pipeline)) +
            geom_col() +
            geom_text(aes(label = round(total_min, 1)),
                      vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
            scale_fill_manual(values = aligner_colours) +
            theme_bw() + theme(legend.position = "none") +
            labs(x = NULL, y = "wall-clock (min)")
        p_m = ggplot(pipe_mem, aes(pipeline, peak_rss_gb, fill = pipeline)) +
            geom_col() +
            geom_text(aes(label = round(peak_rss_gb, 1)),
                      vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
            scale_fill_manual(values = aligner_colours) +
            theme_bw() + theme(legend.position = "none") +
            labs(x = NULL, y = "peak RSS (GB)")
        pp_save_pdf(p_t + p_m, pdir, "bio_perf", width = 6, height = 3)
    }

    panel_theme2 = paper_theme(10)
    panels2 = list()
    if (exists("bc_lists", inherits = FALSE)) {
        ## Smaller internal render dims → UpSetR text uses a larger fraction
        ## of the canvas, so fonts look bigger after patchwork scales the
        ## raster up to fill the grid cell. Title names the sample so
        ## the reader sees at a glance whether this is fig 2 (HeLa) or
        ## fig 3 (mouse skin).
        upset_title = if (is_hela) "Cell-barcode overlap (HeLa)"
                      else         "Cell-barcode overlap (mouse skin)"
        panels2$A = upset_panel(bc_lists, width_in = 5, height_in = 3.5,
                                title = upset_title)
    }
    if (!is.null(per_cell_qc)) {
        qc_long = melt(per_cell_qc, id.vars = "pipeline",
                       measure.vars = c("total_umi", "n_genes", "mito_pct"),
                       variable.name = "metric", value.name = "value")
        metric_labels = c(total_umi = "UMIs\n(log10)",
                          n_genes = "Features\n(log10)",
                          mito_pct = "Mitoc.\n(%)")
        qc_long[, metric := factor(metric, levels = names(metric_labels),
                                   labels = metric_labels)]
        qc_long[metric %in% metric_labels[c("total_umi", "n_genes")],
                value := log10(pmax(value, 1))]
        qc_long[, pipeline := order_pipelines(pipeline, aligners)]
        panels2$B = ggplot(qc_long, aes(pipeline, value, fill = pipeline)) +
            geom_violin(scale = "width", width = 0.8, linewidth = 0.2) +
            geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                         linewidth = 0.3) +
            facet_wrap(~ metric, scales = "free_y", nrow = 1) +
            scale_fill_manual(values = aligner_colours) +
            panel_theme2 + theme(legend.position = "none",
                                 axis.text.x = element_text(angle = 30,
                                                            hjust = 1),
                                 strip.text = element_text(size = 9)) +
            labs(x = "pipeline", y = "per-cell value")
    }
    ## Shared theme tweak for the three square heatmaps (MARD, cluster ARI,
    ## phase ARI): pull the legend in close to the tile grid and keep
    ## panel titles plain (not bold) so they read like a caption.
    heatmap_legend_theme = theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        aspect.ratio = 1,
        plot.title = element_text(size = 9, face = "plain",
                                  margin = margin(b = 2)),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.key.size = grid::unit(0.35, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))
    if (exists("p", inherits = FALSE) &&
        file.exists(biords("pseudobulk_mard"))) {
        panels2$C = p + panel_theme2 + heatmap_legend_theme
    }
    ## Per-cell cross-aligner r density → HeLa panel D.
    if (is_hela && !is.null(p_per_cell_cor)) {
        panels2$D = p_per_cell_cor + panel_theme2 +
            theme(panel.grid = element_blank(),
                  aspect.ratio = 1,
                  plot.title = element_text(size = 10, face = "plain",
                                            margin = margin(b = 2)),
                  plot.subtitle = element_text(size = 8, face = "plain"),
                  legend.position = "right",
                  legend.box.spacing = grid::unit(2, "pt"),
                  legend.key.size = grid::unit(0.4, "cm"))
    }
    ## Cluster ARI heatmap (hela) → panel E.
    if (is_hela && exists("p_ari_heat", inherits = FALSE)) {
        panels2$E = p_ari_heat + panel_theme2 + heatmap_legend_theme
    }
    ## Phase ARI heatmap (hela) → panel F.
    if (is_hela && exists("p_phase_ari", inherits = FALSE)) {
        panels2$F = p_phase_ari + panel_theme2 + heatmap_legend_theme
    }
    ## Cluster UMAP row (hela) → panel G. Wrapped as one patchwork element
    ## so it receives a single tag letter.
    if (is_hela && exists("p_cluster_row", inherits = FALSE)) {
        panels2$G = patchwork::wrap_elements(
            full = patchwork::patchworkGrob(p_cluster_row))
    }
    ## Phase UMAP row (hela) → panel H.
    if (is_hela && exists("p_phase_row", inherits = FALSE)) {
        panels2$H = patchwork::wrap_elements(
            full = patchwork::patchworkGrob(p_phase_row))
    }
    ## Perf: two standalone ggplots so fig 2 gets separate tags I (time)
    ## and J (memory) stacked on the right of the UMAP block. Tight
    ## horizontal margins so the y-axis label sits next to the plot.
    if (!is.null(pipe_time) && !is.null(pipe_mem)) {
        bar_theme2 = panel_theme2 +
            theme(legend.position = "none",
                  plot.margin = grid::unit(c(2, 2, 2, 2), "pt"),
                  axis.title.y = element_text(size = 9,
                                              margin = margin(r = 2)),
                  axis.text.x = element_text(angle = 30, hjust = 1))
        panels2$I = ggplot(pipe_time, aes(pipeline, total_min,
                                          fill = pipeline)) +
            geom_col(width = 0.55) +
            geom_text(aes(label = round(total_min, 1)),
                      vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
            scale_fill_manual(values = aligner_colours) +
            bar_theme2 + labs(x = "pipeline", y = "wall-clock (min)")
        panels2$J = ggplot(pipe_mem, aes(pipeline, peak_rss_gb,
                                         fill = pipeline)) +
            geom_col(width = 0.55) +
            geom_text(aes(label = round(peak_rss_gb, 1)),
                      vjust = -0.3, size = 3) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
            scale_fill_manual(values = aligner_colours) +
            bar_theme2 + labs(x = "pipeline", y = "peak RSS (GB)")
    }
    ## Sendoel-only compositions (legacy; preserved so Fig 3 still builds).
    ## Sendoel UMAPs (celltype and cluster) come from the !is_hela branch
    ## above where p_umap (celltype) and p_cl (cluster) are built.
    if (!is_hela && exists("p_umap", inherits = FALSE)) {
        panels2$sendoel_umap = patchwork::wrap_elements(
            full = patchwork::patchworkGrob(p_umap))
    }
    if (!is_hela && exists("p_cl", inherits = FALSE)) {
        panels2$sendoel_cl = patchwork::wrap_elements(
            full = patchwork::patchworkGrob(p_cl))
    }
    if (!is_hela && exists("p_bcor", inherits = FALSE)) {
        panels2$sendoel_bcor = p_bcor + panel_theme2 +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  aspect.ratio = 1)
    }
    if (!is_hela && exists("p_ari_combo", inherits = FALSE) &&
        !is.null(p_ari_combo)) {
        panels2$sendoel_ari = p_ari_combo + panel_theme2 +
            theme(axis.text.x = element_text(angle = 30, hjust = 1),
                  legend.position = "top")
    }
    blank2 = function() ggplot() + theme_void()
    required_panels = if (is_hela) c("A", "B", "C", "D", "E", "F",
                                     "G", "H", "I", "J")
                      else         c("A", "B", "C")
    for (k in required_panels) {
        if (is.null(panels2[[k]])) {
            warning("fig2 panel ", k, " missing; rendering blank. wd=", wd)
            panels2[[k]] = blank2()
        }
    }
    if (is_hela) {
        ## HeLa fig 2: 18 cols x 8 rows, A-J. Orphan rows removed;
        ## every panel now spans exactly 2 rows so no blank strips
        ## remain under A, G, or H.
        ##   Row 1-2: A UpSet (6), B QC violins (6), C MARD heatmap (6)
        ##   Row 3-4: D per-cell r (6), E cluster ARI (6), F phase ARI (6)
        ##   Row 5-6: G cluster UMAPs (14), I perf wall-clock (4)
        ##   Row 7-8: H phase UMAPs (14), J perf peak RSS (4)
        fig2_design = paste(
            "AAAAAABBBBBBCCCCCC",
            "AAAAAABBBBBBCCCCCC",
            "DDDDDDEEEEEEFFFFFF",
            "DDDDDDEEEEEEFFFFFF",
            "GGGGGGGGGGGGGII###",
            "GGGGGGGGGGGGGII###",
            "HHHHHHHHHHHHHJJ###",
            "HHHHHHHHHHHHHJJ###", sep = "\n")
        fig2 = patchwork::wrap_plots(
            A = panels2$A, B = panels2$B, C = panels2$C,
            D = panels2$D, E = panels2$E, F = panels2$F,
            G = panels2$G, H = panels2$H,
            I = panels2$I, J = panels2$J,
            design = fig2_design,
            heights = c(1.0, 1.0, 1.0, 1.0, 0.97, 0.97, 0.97, 0.97)) +
            patchwork::plot_annotation(tag_levels = "A") &
            paper_shared_theme
        pp_save_pdf(fig2, pdir, fig_stem, width = 18, height = 10)
    } else {
        ## Sendoel fig 3: narrative-ordered A-G. See the is_hela branch above
        ## for the canonical "2/3 UMAPs left, 1/3 perf right" layout; this
        ## branch uses an earlier stacked layout and is scheduled to move to
        ## the same 12-col design as HeLa on the next pass.
        ## Panel I/J are the separate perf ggplots shared with HeLa.
        bar_theme_s = paper_theme(10) +
            theme(legend.position = "none",
                  plot.margin = grid::unit(c(2, 8, 2, 8), "pt"),
                  axis.text.x = element_text(angle = 30, hjust = 1))
        sendoel_perf = (panels2$I + panels2$J) +
            patchwork::plot_layout(ncol = 2)
        fig3_design = paste(
            "AABBBBCCCC",
            "AABBBBCCCC",
            "DDDDDDDDDD",
            "DDDDDDDDDD",
            "EEEEEEEEEE",
            "EEEEEEEEEE",
            "FFFFFGGGGG",
            "FFFFFGGGGG", sep = "\n")
        fig3 = patchwork::wrap_plots(
            A = panels2$A, B = panels2$B, C = panels2$C,
            D = panels2$sendoel_cl,
            E = panels2$sendoel_umap,
            F = panels2$sendoel_ari,
            G = sendoel_perf,
            design = fig3_design,
            heights = c(0.8, 0.8, 1, 1, 1, 1, 0.8, 0.8)) +
            patchwork::plot_annotation(tag_levels = "A") &
            paper_shared_theme
        pp_save_pdf(fig3, pdir, fig_stem, width = 10.5, height = 14)
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

    pipe_time = bm[pipeline %in% aligners & !is_install_rule(file),
                   .(total_min = sum(minutes)), by = pipeline]
    pipe_mem = bm[pipeline %in% aligners & !is_install_rule(file),
                  .(peak_rss_gb = max(max_rss_gb, na.rm = TRUE)), by = pipeline]
    pp_save_csv(pipe_time, pdir, "bench_total_time")
    pp_save_csv(pipe_mem, pdir, "bench_peak_memory")

    p_t = ggplot(pipe_time, aes(pipeline, total_min, fill = pipeline)) +
        geom_col() + geom_text(aes(label = round(total_min, 1)),
                               vjust = -0.3, size = 3) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
        scale_fill_manual(values = aligner_colours) +
        theme_bw() + theme(legend.position = "none") +
        labs(x = NULL, y = "total time (min)")
    p_m = ggplot(pipe_mem, aes(pipeline, peak_rss_gb, fill = pipeline)) +
        geom_col() + geom_text(aes(label = round(peak_rss_gb, 1)),
                               vjust = -0.3, size = 3) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
        scale_fill_manual(values = aligner_colours) +
        theme_bw() + theme(legend.position = "none") +
        labs(x = NULL, y = "peak RSS (GB)")
    pp_save_pdf(p_t + p_m, pdir, paste0(opt$bench_prefix, "_bench_total"),
                width = 6, height = 3)

    bm_steps = bm[pipeline %in% aligners]
    if (nrow(bm_steps) > 0) {
        p_step = ggplot(bm_steps, aes(reorder(step, minutes), minutes,
                                      fill = pipeline)) +
            geom_col(width = 0.7) + coord_flip() +
            scale_fill_manual(values = aligner_colours) +
            theme_bw() +
            labs(x = NULL, y = "wall-clock time (min)", fill = NULL)
        pp_save_pdf(p_step, pdir, paste0(opt$bench_prefix, "_bench_step"),
                    width = 6, height = 4.5)
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

## One-off source-compilation rules matched by benchmark file name. Excluded
## from per-pipeline totals so figures reflect per-run cost rather than the
## first-invocation install cost.
is_install_rule = function(file) grepl("_install\\.txt$", file)

load_pipeline_time = function(bench_dir, aligners) {
    bm = load_benchmarks(bench_dir, aligners)
    bm[pipeline %in% aligners & !is_install_rule(file),
       .(total_min = sum(minutes)), by = pipeline]
}

load_pipeline_memory = function(bench_dir, aligners) {
    bm = load_benchmarks(bench_dir, aligners)
    bm[pipeline %in% aligners & !is_install_rule(file),
       .(peak_rss_gb = max(max_rss_gb, na.rm = TRUE)), by = pipeline]
}

## Per-sample linker error rate panel. Consumes the TSV written by the
## rhapsodist linker_qc rule (scan_r1_linkers over the first 10,000 R1 reads)
## and produces a one-panel PDF showing the Hamming-distance distribution
## against both v1 and enhanced linker templates.
run_linker_qc = function(opt) {
    wd = opt$working_dir
    samp = opt$sample
    tsv = file.path(wd, "linker_qc", paste0(samp, "_linker_errors.tsv"))
    stopifnot(file.exists(tsv))
    pdir = setup_paper_dir(wd, samp)

    raw = readLines(tsv)
    body = raw[!grepl("^#", raw)]
    dt = fread(text = paste(body, collapse = "\n"))
    if (nrow(dt) == 0) {
        warning("empty linker QC TSV, skipping: ", tsv)
        return(invisible(NULL))
    }

    detected = NA_character_
    hit = grep("^# detected_class\\t", raw, value = TRUE)
    if (length(hit) > 0) detected = strsplit(hit[1], "\t")[[1]][2]
    subtitle = sprintf("Auto-detected chemistry: %s", detected)

    max_err = min(max(dt$n_errors), 10)
    pd = dt[n_errors <= max_err]
    p = ggplot(pd, aes(x = factor(n_errors), y = frac, fill = class)) +
        geom_col(position = position_dodge(width = 0.9)) +
        geom_text(aes(label = scales::percent(frac, accuracy = 0.1)),
                  position = position_dodge(width = 0.9),
                  vjust = -0.3, size = 2.4) +
        scale_fill_brewer(palette = "Set2") +
        scale_y_continuous(labels = scales::percent,
                           expand = expansion(mult = c(0.02, 0.12))) +
        theme_bw(base_size = 10) +
        theme(panel.grid = element_blank()) +
        labs(title = paste0("Linker mismatches per read (", samp, ")"),
             subtitle = subtitle,
             x = "Linker mismatches per read (sum across both linkers)",
             y = "Fraction of reads",
             fill = "Candidate class")

    pp_save_pdf(p, pdir, "bio_linker_errors", width = 7, height = 4.5)
    pp_save_csv(dt, pdir, "bio_linker_errors")
    message("linker QC figure written to ", pdir)
}

main = function() {
    opt = parse_args()
    switch(opt$kind,
        simulation = run_simulation(opt),
        biology    = run_biology(opt),
        benchmarks = run_benchmarks(opt),
        linker_qc  = run_linker_qc(opt))
}

if (!interactive()) main()
