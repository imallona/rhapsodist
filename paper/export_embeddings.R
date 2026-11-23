#!/usr/bin/env Rscript
## Writes UMAP coordinates and per-cell labels (cluster, cell type, cell cycle
## phase) of every aligner of one sample as a CSV. Reads the Seurat objects and
## the clusters table saved by 04_biology.Rmd.
##
## Usage:
##   Rscript paper/export_embeddings.R --working_dir output/sendoel2024 \
##       --sample sample_16_wta_p60 --aligners starsolo,kallisto,alevin,sbg \
##       --out bio_umap.csv

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(Seurat)
})

p = ArgumentParser()
p$add_argument("--working_dir", required = TRUE)
p$add_argument("--sample", required = TRUE)
p$add_argument("--aligners", default = "starsolo,kallisto,alevin,sbg")
p$add_argument("--out", required = TRUE)
opt = p$parse_args()

biology_rds = function(stem) {
    file.path(opt$working_dir, sprintf("%s_biology_%s.rds", opt$sample, stem))
}

aligner_umap = function(aligner) {
    so = readRDS(biology_rds(paste0(aligner, "_seurat")))
    emb = Embeddings(so, "umap")
    data.table(pipeline = aligner, barcode = rownames(emb),
               umap1 = emb[, 1], umap2 = emb[, 2])
}

aligners = strsplit(opt$aligners, ",")[[1]]
umap = rbindlist(lapply(aligners, aligner_umap))
labels = as.data.table(readRDS(biology_rds("clusters")))
keep = intersect(c("pipeline", "barcode", "cluster_prefixed", "celltype", "Phase"),
                 colnames(labels))
umap = merge(umap, labels[, ..keep], by = c("pipeline", "barcode"), all.x = TRUE)
fwrite(umap, opt$out)
message(sprintf("%d cells written", nrow(umap)))
