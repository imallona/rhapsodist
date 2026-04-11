## GTF parsing helpers shared by SCE generators.

suppressPackageStartupMessages({
    library(S4Vectors)
})

## Parse a Gencode/Ensembl GTF and return one row per gene feature with
## id (gene_id, as written in the GTF including any version suffix),
## name (gene_name if present, else id) and type (gene_type/gene_biotype if present).
parse_gtf_genes <- function(gtf_file) {
    con <- file(gtf_file, open = "r")
    on.exit(close(con), add = TRUE)

    chunks <- list()
    chunk_size <- 200000L
    repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (length(lines) == 0L) break
        lines <- lines[!startsWith(lines, "#")]
        if (length(lines) == 0L) next
        keep <- grepl("\tgene\t", lines, fixed = TRUE)
        if (any(keep)) chunks[[length(chunks) + 1L]] <- lines[keep]
    }
    if (length(chunks) == 0L) {
        return(S4Vectors::DataFrame(id = character(0),
                                    name = character(0),
                                    type = character(0)))
    }

    gene_lines <- unlist(chunks, use.names = FALSE)
    fields <- strsplit(gene_lines, "\t", fixed = TRUE)
    attrs <- vapply(fields, `[[`, character(1), 9L)

    gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", attrs)

    has_name <- grepl('gene_name "', attrs, fixed = TRUE)
    gene_name <- gene_id
    gene_name[has_name] <- sub('.*gene_name "([^"]+)".*', "\\1", attrs[has_name])

    has_type_gencode <- grepl('gene_type "', attrs, fixed = TRUE)
    has_type_ensembl <- grepl('gene_biotype "', attrs, fixed = TRUE)
    gene_type <- rep(NA_character_, length(attrs))
    gene_type[has_type_gencode] <- sub('.*gene_type "([^"]+)".*', "\\1",
                                       attrs[has_type_gencode])
    gene_type[has_type_ensembl] <- sub('.*gene_biotype "([^"]+)".*', "\\1",
                                       attrs[has_type_ensembl])

    S4Vectors::DataFrame(id = gene_id, name = gene_name, type = gene_type)
}

## Build a rowData DataFrame aligned with the given gene_ids.
## If a gene_id is not found in the GTF table, name falls back to the id itself
## and type is NA. Matching first tries the id as written, then falls back to the
## id with any trailing version suffix stripped (".N"), to handle mismatches
## between versioned and unversioned identifiers.
build_rowdata_from_gtf <- function(gene_ids, gtf_table) {
    idx <- match(gene_ids, gtf_table$id)
    missing_idx <- is.na(idx)
    if (any(missing_idx)) {
        stripped <- sub("\\.[0-9]+$", "", gene_ids[missing_idx])
        gtf_stripped <- sub("\\.[0-9]+$", "", gtf_table$id)
        idx[missing_idx] <- match(stripped, gtf_stripped)
    }

    name <- ifelse(is.na(idx), gene_ids, gtf_table$name[idx])
    type <- ifelse(is.na(idx), NA_character_, gtf_table$type[idx])

    S4Vectors::DataFrame(name = name,
                         type = type,
                         value = NA_character_,
                         row.names = gene_ids)
}
