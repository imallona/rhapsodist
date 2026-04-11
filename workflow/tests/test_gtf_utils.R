#!/usr/bin/env Rscript
## Self-contained tests for workflow/src/gtf_utils.R.
## Run with:
##   cd workflow && Rscript tests/test_gtf_utils.R
## Any failure aborts via stopifnot with a descriptive message.

script_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    f_flag <- which(ca == "-f")
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else if (length(f_flag) > 0 && length(ca) > f_flag[1]) {
        dirname(normalizePath(ca[f_flag[1] + 1L]))
    } else {
        "."
    }
})
source(file.path(script_dir, "..", "src", "gtf_utils.R"))

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

fixture_lines <- c(
    '##description: mini GTF fixture',
    '#!genome-build GRCm39',
    'chr1\tHAVANA\tgene\t3073253\t3074322\t.\t+\t.\tgene_id "ENSMUSG00000102693.2"; gene_type "TEC"; gene_name "4933401J01Rik"; level 2;',
    'chr1\tHAVANA\ttranscript\t3073253\t3074322\t.\t+\t.\tgene_id "ENSMUSG00000102693.2"; transcript_id "ENSMUST00000193812.2"; gene_name "4933401J01Rik";',
    'chr1\tENSEMBL\tgene\t3102016\t3102125\t.\t+\t.\tgene_id "ENSMUSG00000064842.3"; gene_type "snRNA"; gene_name "Gm26206"; level 3;',
    'chrM\tENSEMBL\tgene\t2751\t3707\t.\t+\t.\tgene_id "ENSMUSG00000064341.1"; gene_type "protein_coding"; gene_name "mt-Nd1"; level 3;',
    'chr1\tHAVANA\tgene\t4000000\t4000500\t.\t-\t.\tgene_id "ENSMUSG00000999999.1"; gene_type "lncRNA"; level 2;',
    'chr1\tENSEMBL\texon\t3102016\t3102125\t.\t+\t.\tgene_id "ENSMUSG00000064842.3"; exon_number 1;'
)

fixture_ensembl_lines <- c(
    '#!genome-build GRCm39',
    'chr1\tensembl\tgene\t3073253\t3074322\t.\t+\t.\tgene_id "ENSMUSG00000102693"; gene_version "2"; gene_name "4933401J01Rik"; gene_biotype "TEC";',
    'chr1\tensembl\tgene\t3102016\t3102125\t.\t+\t.\tgene_id "ENSMUSG00000064842"; gene_version "3"; gene_name "Gm26206"; gene_biotype "snRNA";'
)

tmp_gencode <- tempfile(fileext = ".gtf")
tmp_ensembl <- tempfile(fileext = ".gtf")
tmp_empty   <- tempfile(fileext = ".gtf")
on.exit(unlink(c(tmp_gencode, tmp_ensembl, tmp_empty)), add = TRUE)

writeLines(fixture_lines, tmp_gencode)
writeLines(fixture_ensembl_lines, tmp_ensembl)
writeLines(c("## empty", "#!only comments"), tmp_empty)

## parse_gtf_genes tests

gencode_df <- parse_gtf_genes(tmp_gencode)
check("gencode: four gene rows parsed",
      nrow(gencode_df) == 4L)
check("gencode: id column contains versioned ids",
      identical(sort(gencode_df$id),
                sort(c("ENSMUSG00000102693.2", "ENSMUSG00000064842.3",
                       "ENSMUSG00000064341.1", "ENSMUSG00000999999.1"))))
check("gencode: name falls back to id when gene_name is absent",
      gencode_df$name[gencode_df$id == "ENSMUSG00000999999.1"] == "ENSMUSG00000999999.1")
check("gencode: mitochondrial symbol extracted",
      gencode_df$name[gencode_df$id == "ENSMUSG00000064341.1"] == "mt-Nd1")
check("gencode: gene_type captured",
      gencode_df$type[gencode_df$id == "ENSMUSG00000064842.3"] == "snRNA")
check("gencode: transcript/exon rows not included",
      all(gencode_df$id != "ENSMUST00000193812.2"))

ensembl_df <- parse_gtf_genes(tmp_ensembl)
check("ensembl: two gene rows parsed",
      nrow(ensembl_df) == 2L)
check("ensembl: gene_biotype captured as type",
      ensembl_df$type[ensembl_df$id == "ENSMUSG00000064842"] == "snRNA")

empty_df <- parse_gtf_genes(tmp_empty)
check("empty gtf returns a zero-row DataFrame",
      nrow(empty_df) == 0L)
check("empty gtf keeps expected column schema",
      all(c("id", "name", "type") %in% colnames(empty_df)))

## build_rowdata_from_gtf tests

query_exact <- c("ENSMUSG00000102693.2", "ENSMUSG00000064341.1",
                 "ENSMUSG00000999999.1", "ENSMUSG00000000001.9")
rd_exact <- build_rowdata_from_gtf(query_exact, gencode_df)
check("rowdata row order matches query order",
      identical(rownames(rd_exact), query_exact))
check("rowdata: exact match yields gene_name",
      rd_exact$name[1] == "4933401J01Rik")
check("rowdata: mito symbol preserved",
      rd_exact$name[2] == "mt-Nd1")
check("rowdata: missing gene_name in GTF falls back to id",
      rd_exact$name[3] == "ENSMUSG00000999999.1")
check("rowdata: gene id absent from GTF falls back to itself",
      rd_exact$name[4] == "ENSMUSG00000000001.9")
check("rowdata: value column is NA",
      all(is.na(rd_exact$value)))
check("rowdata: columns are name, type, value",
      identical(colnames(rd_exact), c("name", "type", "value")))

query_versioned <- c("ENSMUSG00000102693.7")
rd_version <- build_rowdata_from_gtf(query_versioned, ensembl_df)
check("rowdata: unversioned GTF matches versioned query via fallback",
      rd_version$name[1] == "4933401J01Rik")

## build_rowdata_from_gtf emits exactly one row per query id, even duplicates
query_dup <- c("ENSMUSG00000064341.1", "ENSMUSG00000064341.1")
rd_dup <- build_rowdata_from_gtf(query_dup, gencode_df)
check("rowdata: duplicate query ids produce two rows",
      nrow(rd_dup) == 2L && all(rd_dup$name == "mt-Nd1"))

## mito detection: after an Rmd-style prefix, "^mt-" hits the right rows
combined <- paste0(rd_exact$name, "__", rownames(rd_exact))
check("combined rowname starts with symbol when symbol is known",
      combined[2] == "mt-Nd1__ENSMUSG00000064341.1")
check("mito regex matches mt- symbols in combined names",
      length(grep("^mt-", combined, ignore.case = TRUE)) == 1L)

cat("\nAll gtf_utils tests passed.\n")
