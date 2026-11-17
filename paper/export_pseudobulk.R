#!/usr/bin/env Rscript
## Writes the per-aligner pseudobulk (summed counts on shared gene symbols)
## of one sample as a CSV, as 04_biology.Rmd computes it, and optionally the
## number of cells called by each aligner.
##
## Usage:
##   Rscript paper/export_pseudobulk.R --working_dir output/gse282765 \
##       --sample colon --aligners starsolo,kallisto,alevin,sbg --out pb.csv \
##       --cells_out cells.csv

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(SingleCellExperiment)
})

p = ArgumentParser()
p$add_argument("--working_dir", required = TRUE)
p$add_argument("--sample", required = TRUE)
p$add_argument("--aligners", default = "starsolo,kallisto,alevin,sbg")
p$add_argument("--out", required = TRUE)
p$add_argument("--cells_out", default = "")
opt = p$parse_args()

aligner_pseudobulk = function(aligner) {
    fn = file.path(opt$working_dir, aligner, opt$sample,
                   paste0(opt$sample, "_", aligner, "_sce.rds"))
    if (!file.exists(fn)) return(NULL)
    sce = readRDS(fn)
    sce = sce[rowSums(counts(sce)) > 0, colSums(counts(sce)) > 0]
    rd = rowData(sce)
    sym = rownames(sce)
    for (col in c("name", "SYMBOL", "gene_name")) {
        if (col %in% colnames(rd)) { sym = as.character(rd[[col]]); break }
    }
    rownames(sce) = make.unique(sym)
    n_cells[[aligner]] <<- ncol(sce)
    rowSums(counts(sce))
}

n_cells = list()

aligners = strsplit(opt$aligners, ",")[[1]]
pb = Filter(Negate(is.null), setNames(lapply(aligners, aligner_pseudobulk), aligners))
shared_genes = Reduce(intersect, lapply(pb, names))
pb_dt = data.table(gene = shared_genes,
                   do.call(cbind, lapply(pb, function(v) v[shared_genes])))
fwrite(pb_dt, opt$out)
if (nzchar(opt$cells_out))
    fwrite(data.table(pipeline = names(n_cells), n_cells = unlist(n_cells)),
           opt$cells_out)
message(sprintf("%d shared genes, %d aligners", length(shared_genes), length(pb)))
