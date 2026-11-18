#!/usr/bin/env Rscript
## Draws the manuscript figures at print size (6.3 in wide, 7 pt text) from the
## CSV files written by assemble_paper_figures.R, export_embeddings.R and
## export_pseudobulk.R.
##
## Usage:
##   Rscript paper/manuscript_figures.R --results <dir> --out <dir>
##
## <dir> holds one folder per run (simul, simul_sparse, moro_mallona2025,
## sendoel2024, sendoel2024_downsampled40, gse282765, gse301173), each with
## the paper/ folder of that run.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(ggrastr)
})

p = ArgumentParser()
p$add_argument("--results", required = TRUE)
p$add_argument("--out", required = TRUE)
p$add_argument("--only", default = "", help = "comma-separated figure names")
opt = p$parse_args()
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

script_dir = tryCatch({
    f = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    dirname(normalizePath(f))
}, error = function(e) getwd())
source(file.path(script_dir, "..", "workflow", "src", "agreement_metrics.R"))

page_width = 6.3
base_size = 7
label_size = base_size * 0.8 / ggplot2::.pt

aligner_order = c("starsolo", "kallisto", "alevin", "sbg")
aligner_labels = c(starsolo = "STARsolo", kallisto = "kallisto", alevin = "alevin",
                   sbg = "sbg", truth = "truth")
aligner_colours = c(starsolo = "#E69F00", kallisto = "#56B4E9", alevin = "#009E73",
                    sbg = "#CC79A7", truth = "#000000")

datasets = list(
    hela = list(run = "moro_mallona2025", sample = "hela_unmod", label = "HeLa"),
    epidermis = list(run = "sendoel2024", sample = "sample_16_wta_p60",
                     label = "Mouse epidermis", markers = "skin_markers.tsv"),
    neutrophils = list(run = "gse301173", sample = "neutrophils",
                       label = "Mouse neutrophils", markers = "neutrophil_markers.tsv"),
    colon = list(run = "gse282765", sample = "colon", label = "Mouse colon"))

theme_paper = function() {
    theme_classic(base_size = base_size) +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = base_size),
              axis.text = element_text(size = base_size - 1, colour = "black"),
              axis.line = element_line(linewidth = 0.3),
              axis.ticks = element_line(linewidth = 0.3),
              legend.key.size = grid::unit(3, "mm"),
              legend.text = element_text(size = base_size - 1),
              legend.title = element_text(size = base_size),
              plot.title = element_text(size = base_size, face = "plain",
                                        hjust = 0.5),
              plot.tag = element_text(size = base_size + 2, face = "bold"),
              plot.margin = margin(3, 3, 3, 3))
}
theme_set(theme_paper())

sample_dir = function(key) {
    d = datasets[[key]]
    file.path(opt$results, d$run, "paper", d$sample)
}
read_sample = function(key, name) fread(file.path(sample_dir(key), name))

save_figure = function(plot, name, height, width = page_width) {
    fn = file.path(opt$out, paste0(name, ".pdf"))
    ggsave(fn, plot, width = width, height = height, units = "in",
           device = cairo_pdf)
    message("wrote ", fn)
}

aligner_factor = function(x, levels = aligner_order) {
    factor(x, levels = intersect(c(levels, "truth"), unique(c(levels, x))),
           labels = aligner_labels[intersect(c(levels, "truth"),
                                             unique(c(levels, x)))])
}

## Light fills keep black labels readable; direction -1 when lower is better.
tile_fill = function(limits, direction = 1) {
    colours = c("#F7FBFF", "#9ECAE1")
    if (direction < 0) colours = c("#FFF5EB", "#FDAE6B")
    scale_fill_gradient(low = colours[1], high = colours[2], limits = limits,
                        oob = scales::squish, guide = "none")
}

## Lower triangle of a pairwise statistic, one tile per aligner pair.
pair_heatmap = function(dt, title, digits = 2, limits = c(0, 1),
                        direction = 1, levels = aligner_order) {
    dt = copy(dt)
    setnames(dt, ncol(dt), "value")
    dt = dt[pipeline1 != pipeline2]
    rank1 = match(dt$pipeline1, levels)
    rank2 = match(dt$pipeline2, levels)
    dt = dt[rank1 > rank2]
    dt[, x := aligner_factor(pipeline2, levels)]
    dt[, y := aligner_factor(pipeline1, levels)]
    ggplot(dt, aes(x, y, fill = value)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = formatC(value, format = "f", digits = digits)),
                  size = label_size) +
        tile_fill(limits, direction) +
        scale_y_discrete(limits = rev) +
        coord_fixed() +
        labs(x = NULL, y = NULL, title = title) +
        theme(axis.line = element_blank(), axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
}

## Matrix of a pairwise statistic computed on the pseudobulk table.
pseudobulk_pairs = function(pb, stat, levels = aligner_order) {
    pipes = intersect(levels, colnames(pb))
    grid = CJ(pipeline1 = pipes, pipeline2 = pipes)
    grid[, value := mapply(function(a, b) stat(pb[[a]], pb[[b]]),
                           pipeline1, pipeline2)]
    grid
}

## Upset plot: intersection sizes on top, membership dots below.
upset_plot = function(sets, n_max = 8) {
    pipes = intersect(aligner_order, names(sets))
    all_cells = unique(unlist(sets))
    member = sapply(pipes, function(s) all_cells %in% sets[[s]])
    key = apply(member, 1, function(r) paste(pipes[r], collapse = "&"))
    sizes = as.data.table(table(key))[order(-N)][seq_len(min(.N, n_max))]
    sizes[, rank := factor(seq_len(.N))]
    dots = sizes[, .(pipeline = pipes, present = pipes %in% strsplit(key, "&")[[1]]),
                 by = .(key, rank)]
    dots[, pipeline := aligner_factor(pipeline)]
    bars = ggplot(sizes, aes(rank, N)) +
        geom_col(width = 0.7, fill = "grey30") +
        geom_text(aes(label = format(N, big.mark = ",")), vjust = -0.3,
                  size = label_size * 0.9) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                           labels = scales::label_comma()) +
        labs(x = NULL, y = "cells") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.line.x = element_blank())
    matrix = ggplot(dots, aes(rank, pipeline)) +
        geom_point(aes(colour = present), size = 1.4) +
        geom_line(data = dots[present == TRUE], aes(group = rank),
                  linewidth = 0.4) +
        scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey85"),
                            guide = "none") +
        scale_y_discrete(limits = rev) +
        labs(x = NULL, y = NULL) +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
              axis.line = element_blank())
    wrap_elements(full = bars / matrix + plot_layout(heights = c(1.6, 1)))
}

qc_violins = function(qc) {
    qc = qc[, .(pipeline, total_umi = as.numeric(total_umi),
                n_genes = as.numeric(n_genes), mito_pct = as.numeric(mito_pct))]
    long = melt(qc, id.vars = "pipeline",
                measure.vars = c("total_umi", "n_genes", "mito_pct"),
                variable.name = "metric", value.name = "value")
    long[metric != "mito_pct", value := log10(pmax(value, 1))]
    long[, metric := factor(metric, c("total_umi", "n_genes", "mito_pct"),
                            c("UMIs (log10)", "genes (log10)", "mitochondrial %"))]
    long[, pipeline := aligner_factor(pipeline)]
    ggplot(long, aes(pipeline, value, fill = pipeline)) +
        geom_violin(scale = "width", width = 0.8, linewidth = 0.2) +
        geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                     linewidth = 0.25) +
        facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = setNames(aligner_colours, aligner_labels[names(aligner_colours)]),
                          guide = "none") +
        labs(x = NULL, y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

per_cell_density = function(cor_dt) {
    cor_dt = copy(cor_dt)
    cor_dt[, pair := gsub(" vs ", "-", gsub("starsolo", "STARsolo", pair))]
    ggplot(cor_dt, aes(correlation, colour = pair)) +
        geom_density(linewidth = 0.4) +
        scale_colour_brewer(palette = "Dark2", name = NULL) +
        labs(x = "per-cell Pearson r", y = "density") +
        theme(legend.position = "bottom",
              legend.key.size = grid::unit(2, "mm"),
              legend.text = element_text(size = base_size - 2),
              legend.margin = margin(0, 0, 0, 0)) +
        guides(colour = guide_legend(ncol = 2))
}

umap_panel = function(umap, colour_by, palette = NULL, legend = TRUE) {
    umap = copy(umap)
    umap[, pipeline := aligner_factor(pipeline)]
    umap[, label := get(colour_by)]
    umap = umap[sample(.N)]
    p = ggplot(umap, aes(umap1, umap2, colour = label)) +
        geom_point_rast(size = 0.25, stroke = 0, raster.dpi = 400) +
        facet_wrap(~ pipeline, nrow = 1, scales = "free") +
        labs(x = "UMAP 1", y = "UMAP 2", colour = NULL) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              legend.position = if (legend) "bottom" else "none") +
        guides(colour = guide_legend(override.aes = list(size = 1.5),
                                     nrow = ceiling(uniqueN(umap$label) / 4)))
    if (!is.null(palette)) p = p + scale_colour_manual(values = palette)
    p
}

## Louvain clusters have per-aligner labels; colour them by rank within aligner.
cluster_rank = function(umap) {
    umap = copy(umap)
    umap[, cluster_rank := frank(-tabulate(factor(cluster_prefixed))[
        as.integer(factor(cluster_prefixed))], ties.method = "dense"),
        by = pipeline]
    umap[, cluster_rank := factor(cluster_rank)]
    umap
}

time_memory_bars = function(time_seg, memory) {
    time_seg = copy(time_seg)
    time_seg[, segment := factor(segment, c("preprocessing", "aligner"))]
    totals = time_seg[, .(minutes = sum(minutes)), by = pipeline]
    time_seg[, pipeline := aligner_factor(pipeline)]
    totals[, pipeline := aligner_factor(pipeline)]
    memory = copy(memory)
    memory[, pipeline := aligner_factor(pipeline)]
    fills = setNames(aligner_colours, aligner_labels[names(aligner_colours)])
    time_seg[, fill := ifelse(segment == "preprocessing", "barcode preprocessing", as.character(pipeline))]
    p_time = ggplot(time_seg, aes(minutes, pipeline, fill = fill)) +
        geom_col(width = 0.65, position = position_stack(reverse = TRUE)) +
        geom_text(data = totals, aes(minutes, pipeline,
                                     label = formatC(minutes, format = "f", digits = 1)),
                  inherit.aes = FALSE, hjust = -0.15, size = label_size) +
        scale_fill_manual(values = c(fills, "barcode preprocessing" = "grey60"),
                          breaks = "barcode preprocessing", name = NULL) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
        scale_y_discrete(limits = rev) +
        labs(x = "wall-clock time (min)", y = NULL) +
        theme(legend.position = "inside", legend.position.inside = c(0.9, 0.6))
    p_mem = ggplot(memory, aes(peak_rss_gb, pipeline, fill = pipeline)) +
        geom_col(width = 0.65) +
        geom_text(aes(label = formatC(peak_rss_gb, format = "f", digits = 1)),
                  hjust = -0.15, size = label_size) +
        scale_fill_manual(values = fills, guide = "none") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
        scale_y_discrete(limits = rev) +
        labs(x = "peak memory (GiB)", y = NULL)
    list(time = p_time, memory = p_mem)
}

## Figure 1, panels B to G: simulation.
figure_simulation = function() {
    sim = file.path(opt$results, "simul", "paper", "simulated")
    bench = file.path(opt$results, "simul", "paper")
    truth = fread(file.path(sim, "sim_ground_truth.csv"))
    truth[, pipeline := aligner_factor(pipeline)]
    p_b = ggplot(truth, aes(pipeline, n_recovered, fill = pipeline)) +
        geom_col(width = 0.65) +
        geom_hline(yintercept = truth$n_truth[1], linetype = "dashed",
                   linewidth = 0.3) +
        scale_fill_manual(values = setNames(aligner_colours, aligner_labels[names(aligner_colours)]),
                          guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        labs(x = NULL, y = "true barcodes recovered",
             title = sprintf("recall %.2f, precision %.2f",
                             min(truth$recall), min(truth$precision))) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    tags = fread(file.path(sim, "sim_sampletag_confusion.csv"))
    tags[, `:=`(truth_tag = sub("mouse_sampletag_", "tag ", truth_tag),
                pred_tag = sub("mouse_sampletag_", "tag ", pred_tag))]
    p_c = ggplot(tags, aes(pred_tag, truth_tag)) +
        geom_tile(colour = "white", fill = "#9ECAE1") +
        geom_text(aes(label = N), size = label_size) +
        scale_y_discrete(limits = rev) +
        coord_fixed() +
        labs(x = "called tag", y = "true tag") +
        theme(axis.line = element_blank(), axis.ticks = element_blank())
    sim_levels = c(aligner_order, "truth")
    mard = fread(file.path(sim, "sim_pseudobulk_mard.csv"))[, .(pipeline1, pipeline2, mard = mard * 100)]
    p_d = pair_heatmap(mard, "pseudobulk MARD (%)", digits = 2,
                       limits = c(0, 1), direction = -1, levels = sim_levels)
    pcor = fread(file.path(sim, "sim_pseudobulk_correlation.csv"))
    p_g = pair_heatmap(pcor, "pseudobulk Pearson r", digits = 3,
                       limits = c(0.99, 1), levels = sim_levels)
    bars = time_memory_bars(fread(file.path(bench, "bench_total_time_segments.csv")),
                            fread(file.path(bench, "bench_peak_memory.csv")))
    design = "BBCCCDDDGGG\nEEEEEEFFFFF"
    fig = wrap_plots(B = p_b, C = p_c, D = p_d, E = bars$time, F = bars$memory,
                     G = p_g, design = design, heights = c(1.4, 0.8)) +
        plot_annotation(tag_levels = list(c("B", "C", "D", "E", "F", "G")))
    save_figure(fig, "fig1_simulations_panels", height = 3.4)
}

## Figures 2 and 3 and the neutrophil figure share one layout.
figure_biology = function(key, name, second_label, second_title) {
    umap = read_sample(key, "bio_umap.csv")
    sets = split(umap$barcode, umap$pipeline)
    qc = read_sample(key, "bio_per_cell_qc.csv")
    mard = read_sample(key, "bio_pseudobulk_mard.csv")[, .(pipeline1, pipeline2, mard = mard * 100)]
    cor_dt = readRDS(file.path(sample_dir(key), "bio_per_cell_correlation.rds"))
    cluster_ari = read_sample(key, "bio_cluster_ari_matrix.csv")
    second_ari = read_sample(key, switch(second_label,
                                         Phase = "bio_phase_ari_matrix.csv",
                                         "bio_celltype_ari_matrix.csv"))
    run_dir = file.path(opt$results, datasets[[key]]$run, "paper")
    bars = time_memory_bars(fread(file.path(run_dir, "bench_total_time_segments.csv")),
                            fread(file.path(run_dir, "bench_peak_memory.csv")))
    umap = cluster_rank(umap)
    n_clusters = nlevels(umap$cluster_rank)
    cluster_palette = setNames(grDevices::hcl.colors(n_clusters, "Dynamic"),
                               levels(umap$cluster_rank))
    second_levels = sort(unique(na.omit(umap[[second_label]])))
    second_palette = setNames(
        c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
          "#CC79A7", "#AA4499", "#882255", "#44AA99", "#117733",
          "#332288")[seq_along(second_levels)], second_levels)
    if ("none" %in% second_levels) second_palette["none"] = "grey80"
    panels = list(
        A = upset_plot(sets),
        B = qc_violins(qc),
        C = pair_heatmap(mard, "pseudobulk MARD (%)", digits = 0,
                         limits = c(0, 100), direction = -1),
        D = per_cell_density(cor_dt),
        E = pair_heatmap(cluster_ari, "cluster ARI"),
        F = pair_heatmap(second_ari, second_title),
        G = umap_panel(umap, "cluster_rank", cluster_palette, legend = FALSE),
        H = umap_panel(umap, second_label, second_palette),
        I = bars$time,
        J = bars$memory)
    design = "AAAAAABBBBBBB\nCCCDDDDEEEFFF\nGGGGGGGGGGGGG\nHHHHHHHHHHHHH\nIIIIIIJJJJJJJ"
    heights = c(1.15, 1.25, 1.05, 1.35, 0.75)
    if (!is.null(extra <- figure_extra(key))) {
        panels$K = extra
        design = "AAAAAABBBBBBB\nCCCDDDDEEEFFF\nGGGGGGGGGGGGG\nHHHHHHHHHHHHH\nIIIIIJJJJJKKK"
        heights[5] = 1.1
    }
    fig = wrap_plots(panels, design = design, heights = heights) +
        plot_annotation(tag_levels = "A")
    save_figure(fig, name, height = if (is.null(panels$K)) 8.4 else 8.7)
}

## Pseudotime agreement, neutrophils only.
figure_extra = function(key) {
    fn = file.path(sample_dir(key), "bio_pseudotime_spearman.csv")
    if (!file.exists(fn)) return(NULL)
    pair_heatmap(fread(fn), "pseudotime\nSpearman")
}

## Linker mismatch distributions of every experimental dataset.
figure_linker_errors = function() {
    rows = rbindlist(lapply(names(datasets), function(key) {
        dt = read_sample(key, "bio_linker_errors.csv")
        dt[, dataset := datasets[[key]]$label]
    }))
    rows[, dataset := factor(dataset, sapply(datasets, `[[`, "label"))]
    rows[, class := factor(class, c("enhanced", "v1"), c("enhanced linkers", "v1 linkers"))]
    rows = rows[n_errors <= 6]
    p = ggplot(rows, aes(factor(n_errors), frac * 100, fill = class)) +
        geom_col(position = position_dodge(width = 0.8, preserve = "single"),
                 width = 0.75) +
        geom_text(aes(label = ifelse(frac >= 0.005, sprintf("%.1f", frac * 100), "")),
                  position = position_dodge(width = 0.8, preserve = "single"),
                  vjust = -0.3, size = label_size * 0.85) +
        facet_wrap(~ dataset, nrow = 2) +
        scale_fill_manual(values = c("#0072B2", "#D55E00"), name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = "linker mismatches (Hamming distance)", y = "R1 reads (%)") +
        theme(legend.position = "top")
    save_figure(p, "linker_errors", height = 3.6)
}

## Sparse simulation: fraction of entries detected, by true count.
figure_sim_sparse = function() {
    dt = fread(file.path(opt$results, "simul_sparse", "paper", "simulated",
                         "sim_detection_by_true_count.csv"))
    dt = dt[n_entries > 0]
    dt[, true_count := factor(true_count, unique(true_count))]
    dt[, pipeline := aligner_factor(pipeline)]
    dt[, label := ifelse(true_count == "0",
                         sprintf("%.2f", detected_fraction * 100),
                         sprintf("%.2f", detected_fraction * 100))]
    entries = unique(dt[, .(true_count, n_entries)])
    p = ggplot(dt, aes(true_count, pipeline, fill = detected_fraction)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = label), size = label_size) +
        tile_fill(c(0, 1)) +
        scale_y_discrete(limits = rev) +
        scale_x_discrete(labels = sprintf("%s\n(n = %s)", entries$true_count,
                                          formatC(entries$n_entries, big.mark = ",",
                                                  format = "d"))) +
        labs(x = "true count of the cell-gene entry", y = NULL,
             title = "entries with a non-zero count (%)") +
        theme(axis.line = element_blank(), axis.ticks = element_blank())
    save_figure(p, "sim_sparse_detection", height = 1.8, width = 5)
}

marker_genes = function(key) {
    fn = datasets[[key]]$markers
    if (is.null(fn)) return(character(0))
    fread(file.path(script_dir, "..", "workflow", "data", "markers", fn))$marker
}

## Pairwise pseudobulk scatterplots, one row per dataset.
figure_pseudobulk_scatter = function() {
    rows = rbindlist(lapply(names(datasets), function(key) {
        pb = read_sample(key, "bio_pseudobulk.csv")
        markers = marker_genes(key)
        pipes = intersect(aligner_order, colnames(pb))
        pairs = combn(pipes, 2, simplify = FALSE)
        rbindlist(lapply(pairs, function(pr) {
            data.table(dataset = datasets[[key]]$label,
                       is_marker = pb$gene %in% markers,
                       pair = sprintf("%s vs\n%s", aligner_labels[pr[1]], aligner_labels[pr[2]]),
                       x = log10(pb[[pr[1]]] + 1), y = log10(pb[[pr[2]]] + 1),
                       r = cor(log10(pb[[pr[1]]] + 1), log10(pb[[pr[2]]] + 1)))
        }))
    }))
    rows[, dataset := factor(dataset, sapply(datasets, `[[`, "label"))]
    rows[, pair := factor(pair, unique(pair))]
    labels = unique(rows[, .(dataset, pair, r)])
    p = ggplot(rows, aes(x, y)) +
        geom_point_rast(data = rows[is_marker == FALSE], size = 0.1, alpha = 0.3,
                        raster.dpi = 400) +
        geom_point(data = rows[is_marker == TRUE], size = 0.4, colour = "#E69F00") +
        geom_abline(linetype = "dashed", linewidth = 0.3, colour = "red") +
        geom_text(data = labels, aes(x = -Inf, y = Inf, label = sprintf("r = %.3f", r)),
                  hjust = -0.1, vjust = 1.4, size = label_size) +
        facet_grid(dataset ~ pair) +
        coord_fixed() +
        labs(x = "log10(pseudobulk count + 1), first aligner",
             y = "log10(pseudobulk count + 1), second aligner") +
        theme(strip.text.y = element_text(angle = 270),
              panel.border = element_rect(fill = NA, linewidth = 0.3))
    save_figure(p, "pseudobulk_scatter", height = 5.2)
    labels
}

## Spearman correlation and scaled MARD on the pseudobulk, per dataset.
figure_scaled_agreement = function() {
    panels = list()
    for (key in names(datasets)) {
        pb = read_sample(key, "bio_pseudobulk.csv")
        label = datasets[[key]]$label
        sp = pseudobulk_pairs(pb, spearman)
        sm = pseudobulk_pairs(pb, scaled_mard)[, value := value * 100]
        fwrite(rbind(sp[, stat := "spearman"], sm[, stat := "scaled_mard_pct"]),
               file.path(opt$out, sprintf("pseudobulk_agreement_%s.csv", key)))
        panels[[length(panels) + 1]] = pair_heatmap(sp[, .(pipeline1, pipeline2, value)],
                                                    paste(label, "Spearman", sep = "\n"),
                                                    limits = c(0.8, 1))
        panels[[length(panels) + 1]] = pair_heatmap(sm[, .(pipeline1, pipeline2, value)],
                                                    paste(label, "scaled MARD (%)", sep = "\n"),
                                                    digits = 0, limits = c(0, 80),
                                                    direction = -1)
    }
    fig = wrap_plots(panels, ncol = 4) + plot_annotation(tag_levels = "A")
    save_figure(fig, "pseudobulk_scaled_agreement", height = 3.6)
}

## Effect of 1 and 2 tolerated linker mismatches, relative to 0.
figure_linker_tolerance = function() {
    rows = rbindlist(lapply(c("hela", "epidermis"), function(key) {
        dt = read_sample(key, "linker_tolerance.csv")
        base = dt[max_errors == 0, .(aligner, cells0 = n_cells, umis0 = median_umis)]
        dt = merge(dt[max_errors > 0], base, by = "aligner")
        dt[, dataset := datasets[[key]]$label]
        dt[, .(dataset, aligner, max_errors,
               `cells called, change (%)` = (n_cells / cells0 - 1) * 100,
               `median UMIs per cell, change (%)` = (median_umis / umis0 - 1) * 100,
               `scaled MARD to 0 mismatches (%)` = scaled_mard * 100)]
    }))
    long = melt(rows, id.vars = c("dataset", "aligner", "max_errors"),
                variable.name = "metric")
    long[, aligner := aligner_factor(aligner)]
    long[, mismatches := factor(max_errors)]
    long[, dataset := factor(dataset, c("HeLa", "Mouse epidermis"))]
    p = ggplot(long, aes(aligner, value, fill = mismatches)) +
        geom_hline(yintercept = 0, linewidth = 0.3) +
        geom_col(position = position_dodge(width = 0.75), width = 0.7) +
        geom_text(aes(label = sprintf("%.1f", value),
                      vjust = ifelse(value >= 0, -0.3, 1.2)),
                  position = position_dodge(width = 0.75),
                  size = label_size * 0.85) +
        facet_grid(dataset ~ metric, scales = "free_y") +
        scale_fill_manual(values = c(`1` = "#9ECAE1", `2` = "#2171B5"),
                          name = "tolerated linker mismatches") +
        scale_y_continuous(expand = expansion(mult = 0.15)) +
        labs(x = NULL, y = NULL) +
        theme(legend.position = "top",
              strip.text.y = element_text(angle = 270))
    save_figure(p, "linker_tolerance", height = 3.6)
}

## Sample tag calls on the colon against the tags deposited by the authors.
figure_colon_sampletags = function() {
    dt = read_sample("colon", "sampletag_concordance.csv")
    dt[, deposited_tag := sub("mouse_sampletag_", "tag ", deposited_tag)]
    dt[, outcome := factor(outcome, c("same tag", "multiplet", "other tag",
                                      "undetermined", "absent"))]
    p = ggplot(dt, aes(outcome, deposited_tag, fill = fraction)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(N, big.mark = ","),
                                      fraction * 100),
                      colour = fraction > 0.5), size = label_size, lineheight = 0.9) +
        scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"),
                            guide = "none") +
        scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                             guide = "none") +
        scale_y_discrete(limits = rev) +
        labs(x = "rhapsodist call", y = "deposited tag") +
        theme(axis.line = element_blank(), axis.ticks = element_blank())
    save_figure(p, "colon_sampletags", height = 1.5, width = 4.5)
}

## rhapsodist against UniverSC, zUMIs and OpenPipelines on the epidermis.
figure_workflow_comparison = function() {
    wd = file.path(opt$results, "sendoel2024_downsampled40", "paper")
    cmp = fread(file.path(wd, "sample_16_wta_p60", "workflow_comparison.csv"))
    others = c(universc = "UniverSC", zumis = "zUMIs", openpipelines = "OpenPipelines")
    rh = Reduce(function(x, y) merge(x, y, by = "pipeline"), list(
        fread(file.path(wd, "bench_total_time.csv")),
        fread(file.path(wd, "bench_peak_memory.csv")),
        fread(file.path(wd, "sample_16_wta_p60", "rhapsodist_cells.csv"))))
    rh = rh[match(aligner_order, pipeline)]
    tools = rbind(
        rh[, .(tool = aligner_labels[pipeline], group = "rhapsodist",
               cells = n_cells, hours = total_min / 60, memory = peak_rss_gb)],
        unique(cmp[, .(tool = others[workflow], group = "other workflow",
                       cells = n_cells, hours = seconds / 3600,
                       memory = max_rss_mb / 1024)]))
    tools[, tool := factor(tool, tool)]
    bar = function(column, ylab, digits) {
        ggplot(tools, aes(tool, .data[[column]], fill = group)) +
            geom_col(width = 0.7) +
            geom_text(aes(label = formatC(.data[[column]], format = "f",
                                          digits = digits, big.mark = ",")),
                      vjust = -0.3, size = label_size) +
            scale_fill_manual(values = c(rhapsodist = "#0072B2",
                                         `other workflow` = "grey55"), name = NULL) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                               labels = scales::label_comma()) +
            labs(x = NULL, y = ylab) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    sp = cmp[, .(workflow = factor(others[workflow], others),
                 aligner = aligner_factor(aligner), spearman)]
    p_sp = ggplot(sp, aes(aligner, workflow, fill = spearman)) +
        geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%.2f", spearman)), size = label_size) +
        tile_fill(c(0.7, 1)) +
        scale_y_discrete(limits = rev) +
        coord_fixed() +
        labs(x = "rhapsodist aligner", y = NULL, title = "pseudobulk Spearman") +
        theme(axis.line = element_blank(), axis.ticks = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    fig = wrap_plots(A = bar("cells", "cells called", 0), B = p_sp,
                     C = bar("hours", "wall-clock time (h)", 1),
                     D = bar("memory", "peak memory (GiB)", 1),
                     design = "AB\nCD", widths = c(1.4, 1)) +
        plot_layout(guides = "collect") +
        plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
    save_figure(fig, "workflow_comparison", height = 4.4)
}

## Pairwise confusion matrices of one label, Hungarian-ordered.
figure_confusion = function(key, column, name, log_fill) {
    umap = read_sample(key, "bio_umap.csv")
    pipes = intersect(aligner_order, unique(umap$pipeline))
    ## Axis titles name the aligner: drop cluster prefixes and cell type notes.
    short_label = function(x) sub(" \\(.*", "", sub("^[a-z]+_", "", x))
    pairs = combn(pipes, 2, simplify = FALSE)
    tables = lapply(pairs, function(pr) {
        a = umap[pipeline == pr[1], .(barcode, la = short_label(get(column)))]
        b = umap[pipeline == pr[2], .(barcode, lb = short_label(get(column)))]
        table(merge(a, b, by = "barcode")[, .(la, lb)])
    })
    max_cells = max(sapply(tables, max))
    panels = Map(function(pr, tab) {
        rows = rownames(tab)
        if (!anyNA(suppressWarnings(as.numeric(rows)))) tab = tab[order(as.numeric(rows)), ]
        n = max(dim(tab))
        padded = matrix(0, n, n)
        padded[seq_len(nrow(tab)), seq_len(ncol(tab))] = tab
        match_cols = as.integer(clue::solve_LSAP(padded, maximum = TRUE))
        col_order = match_cols[match_cols <= ncol(tab)]
        col_order = c(col_order, setdiff(seq_len(ncol(tab)), col_order))
        dt = as.data.table(tab)
        dt[, la := factor(la, rownames(tab))]
        dt[, lb := factor(lb, colnames(tab)[col_order])]
        ggplot(dt, aes(lb, la, fill = N)) +
            geom_tile() +
            scale_fill_distiller(palette = "Blues", direction = 1, name = "cells",
                                 trans = if (log_fill) "log1p" else "identity",
                                 limits = c(0, max_cells), labels = scales::comma,
                                 breaks = if (log_fill) c(0, 10, 100, 1000, 10000) else waiver()) +
            scale_y_discrete(limits = rev) +
            labs(x = aligner_labels[pr[2]], y = aligner_labels[pr[1]]) +
            theme(axis.text = element_text(size = base_size - 1),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  axis.line = element_blank(), axis.ticks = element_blank(),
                  aspect.ratio = 1)
    }, pairs, tables)
    fig = wrap_plots(panels, ncol = 3, guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"))
    save_figure(fig, name, height = 5.2)
}

figures = list(
    fig1 = function() figure_simulation(),
    fig2 = function() figure_biology("hela", "fig2_hela", "Phase", "phase ARI"),
    fig3 = function() figure_biology("epidermis", "fig3_epidermis", "celltype", "cell type ARI"),
    neutrophils = function() figure_biology("neutrophils", "neutrophils", "celltype",
                                            "cell type ARI"),
    linker_errors = figure_linker_errors,
    sim_sparse = figure_sim_sparse,
    pseudobulk_scatter = figure_pseudobulk_scatter,
    scaled_agreement = figure_scaled_agreement,
    linker_tolerance = figure_linker_tolerance,
    colon_sampletags = figure_colon_sampletags,
    workflow_comparison = figure_workflow_comparison,
    cluster_confusion = function() figure_confusion("epidermis", "cluster_prefixed",
                                                    "cluster_confusion", log_fill = FALSE),
    celltype_confusion = function() figure_confusion("epidermis", "celltype",
                                                     "celltype_confusion", log_fill = TRUE))
selected = if (nzchar(opt$only)) strsplit(opt$only, ",")[[1]] else names(figures)
for (fig in selected) figures[[fig]]()
