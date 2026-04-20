# Changelog

## Unreleased

- Sample config vocabulary: replaced `bead_version` with `allowedlist` (96/384) and `diversity_insets` (yes/no). Both are optional; bead chemistry is auto-detected from R1 linkers and the declared fields are used as a QC check and for logging.
- Added per-sample `use_sampletags` (yes/no). Sampletag demultiplexing is gated per-sample; `species` is only required when `use_sampletags` is yes.
- Added `cb_umi_max_errors` config key (integer, default 0) exposing cutadapt `-e` for the R1 linker trim.
- Added optional paired-fastq `downsample` (percentage in (0, 100], default 100) with `downsample_seed`. Uses seqtk; per-sample `uses.downsample` overrides the global value.
- Added per-sample linker QC report (`{sample}_linker_qc.html` under `linker_qc/`). The existing bead-class scan is reused to emit a per-read hamming-distance histogram for both v1 and enhanced chemistries. Helps pick a value for `cb_umi_max_errors`.
- Added cross-aligner per-cell correlation distribution plot in the per-sample comparison report (Pearson on log1p counts, one density per aligner pair, matched by cell barcode).

## v0.1.0

First tagged release of rhapsodist, a Snakemake workflow for single-cell BD Rhapsody data. 

It runs starsolo, kallisto and alevin-fry on the same inputs (and optionally the BD Rhapsody pipeline in a Singularity container) from either SRA accessions or local fastq paths, standardises barcodes across aligners, and produces per-sample reports covering QC, pseudobulk agreement (MARD and correlation), per-barcode UMI concordance (Bland-Altman), barcode overlap, per-pipeline Seurat clustering and UMAPs, marker-based cell-type assignment, marker co-expression (UpSet and combination bars), marker dotplots, cluster and cell-type adjusted Rand indices, guide RNA join rates for pooled CRISPR screens, and wall-clock / peak-memory benchmarks per pipeline step. Includes a simulated dataset for CI and integration tests.

Rhapsodist can be used to benchmark/compare aligners/pipelines, or to run one of them.
