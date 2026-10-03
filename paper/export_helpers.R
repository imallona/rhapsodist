## Helpers used by assemble_paper_figures.R to write manuscript-ready CSV and
## PDF artefacts under <working_dir>/paper/<sample>/. Point-heavy scatter
## layers are rasterised via ggrastr so PDF files stay small while axes and
## text remain vector.

setup_paper_dir = function(working_dir, sample = NULL) {
    sub = if (is.null(sample) || !nzchar(sample)) "paper" else file.path("paper", sample)
    pdir = file.path(working_dir, sub)
    dir.create(pdir, showWarnings = FALSE, recursive = TRUE)
    pdir
}

pp_save_csv = function(df, pdir, name) {
    if (is.null(df)) return(invisible(NULL))
    fn = file.path(pdir, paste0(name, ".csv"))
    utils::write.csv(df, fn, row.names = FALSE)
    invisible(fn)
}

pp_save_pdf = function(plot, pdir, name, width = 5, height = 4) {
    if (is.null(plot)) return(invisible(NULL))
    fn = file.path(pdir, paste0(name, ".pdf"))
    ggplot2::ggsave(fn, plot = plot, width = width, height = height,
                    device = grDevices::cairo_pdf)
    invisible(fn)
}

pp_save_base_pdf = function(expr, pdir, name, width = 7, height = 5,
                            keep_last_page = FALSE) {
    caller_env = parent.frame()
    fn = file.path(pdir, paste0(name, ".pdf"))
    grDevices::cairo_pdf(fn, width = width, height = height)
    tryCatch(eval(expr, envir = caller_env), finally = grDevices::dev.off())
    ## UpSetR calls grid.newpage() internally, which leaves a blank page 1 on
    ## cairo_pdf. Keep only page 2 using Ghostscript when available.
    if (isTRUE(keep_last_page) && nzchar(Sys.which("gs"))) {
        tmp = paste0(fn, ".tmp")
        status = suppressWarnings(system2("gs",
            c("-sDEVICE=pdfwrite", "-dNOPAUSE", "-dBATCH", "-dQUIET",
              "-dFirstPage=2", "-dLastPage=2",
              paste0("-sOutputFile=", tmp), fn),
            stdout = FALSE, stderr = FALSE))
        if (status == 0 && file.exists(tmp) && file.info(tmp)$size > 0) {
            file.rename(tmp, fn)
        } else if (file.exists(tmp)) {
            file.remove(tmp)
        }
    }
    invisible(fn)
}

pp_save_latex_table = function(df, pdir, name, caption = NULL, label = NULL,
                               digits = 3) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
    fn = file.path(pdir, paste0(name, ".tex"))
    tex = knitr::kable(df, format = "latex", booktabs = TRUE,
                       digits = digits, caption = caption, label = label)
    writeLines(as.character(tex), fn)
    invisible(fn)
}

pp_rasterise = function(layer, dpi = 300) {
    if (requireNamespace("ggrastr", quietly = TRUE)) {
        ggrastr::rasterise(layer, dpi = dpi)
    } else {
        layer
    }
}

pp_bootstrap_ci = function(x, stat = mean, n_boot = 1000,
                           probs = c(0.025, 0.975), seed = 1) {
    if (length(x) == 0) return(c(NA_real_, NA_real_))
    set.seed(seed)
    b = vapply(seq_len(n_boot),
               function(i) stat(sample(x, replace = TRUE)),
               numeric(1))
    unname(stats::quantile(b, probs = probs, na.rm = TRUE))
}

pp_bootstrap_pair = function(values_a, values_b, stat, n_boot = 1000,
                             probs = c(0.025, 0.975), seed = 1) {
    if (length(values_a) != length(values_b) || length(values_a) == 0) {
        return(c(NA_real_, NA_real_))
    }
    set.seed(seed)
    n = length(values_a)
    b = vapply(seq_len(n_boot), function(i) {
        idx = sample.int(n, replace = TRUE)
        stat(values_a[idx], values_b[idx])
    }, numeric(1))
    unname(stats::quantile(b, probs = probs, na.rm = TRUE))
}

## One row per pair of columns of pb_mat: the statistic and its bootstrap 95%
## CI over genes; value is named by value_name.
pp_pairwise_table = function(pb_mat, stat, value_name, n_boot = 1000) {
    pairs = expand.grid(a = colnames(pb_mat), b = colnames(pb_mat),
                        stringsAsFactors = FALSE)
    out = rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
        a = pairs$a[i]; b = pairs$b[i]
        ci = if (a == b) c(NA, NA) else
            pp_bootstrap_pair(pb_mat[, a], pb_mat[, b], stat, n_boot = n_boot)
        data.table(pipeline1 = a, pipeline2 = b, value = stat(pb_mat[, a], pb_mat[, b]),
                   ci_lo = ci[1], ci_hi = ci[2])
    }))
    setnames(out, "value", value_name)
    out
}

## Heatmap of a pp_pairwise_table column.
pp_pairwise_heatmap = function(dt, value_name, title, label_format, limits) {
    ggplot(dt, aes(pipeline1, pipeline2, fill = .data[[value_name]])) +
        geom_tile(colour = "white") +
        geom_text(aes(label = sprintf(label_format, .data[[value_name]])), size = 2.8) +
        scale_fill_viridis_c(option = "viridis", limits = limits,
                             oob = scales::squish, alpha = 0.75) +
        scale_y_discrete(limits = rev) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1), aspect.ratio = 1) +
        labs(x = NULL, y = NULL, fill = value_name, title = title)
}

## Scaled MARD and Spearman tables and heatmaps, saved as <prefix>_scaled_mard
## and <prefix>_spearman.
pp_save_rank_agreement = function(pb_mat, pdir, prefix) {
    scaled = pp_pairwise_table(pb_mat, scaled_mard, "scaled_mard")
    pp_save_csv(scaled, pdir, paste0(prefix, "_scaled_mard"))
    pp_save_pdf(pp_pairwise_heatmap(scaled, "scaled_mard", "Pseudobulk MARD on counts per million",
                                    "%.3f", c(0, 1)),
                pdir, paste0(prefix, "_scaled_mard"), width = 4.5, height = 3.5)
    rho = pp_pairwise_table(pb_mat, spearman, "spearman")
    pp_save_csv(rho, pdir, paste0(prefix, "_spearman"))
    pp_save_pdf(pp_pairwise_heatmap(rho, "spearman", "Pseudobulk Spearman rho", "%.3f", c(NA, 1)),
                pdir, paste0(prefix, "_spearman"), width = 4.5, height = 3.5)
}
