# Changelog

## Unreleased

- `get_txp2gene` now reads the `transcript_id` and `gene_id` GTF attributes by name with gawk rather than by fixed column position, so GTFs with a different attribute order work. `gtf_origin` no longer affects this rule; it still sets the transcriptome fasta header convention. Added gawk to the salmon conda env.
- Added `alevin_usa` config key. When true, a spliced+unspliced (spliceu) reference is built with pyroe from the genome and GTF and alevin-fry quantifies in USA mode (triggered by the 3-column t2g). The alevin SingleCellExperiment keeps spliced plus ambiguous as the main `counts` assay and adds `spliced`, `unspliced` and `ambiguous` assays. Requires `alevin_sketch: true`.
- Added `cell_filtering: none` to keep every observed barcode (no cell filter) across STARsolo, alevin and kallisto. For STARsolo it overrides `soloCellFilter` to None; for alevin and kallisto the knee filter is skipped.
- The transcriptome input now accepts a plain `.fa` as well as `.fa.gz` (deversion uses `gzip -dcf`). The previous `zcat` failed on uncompressed fasta.
- Documented that BD sample tags are extracted from the WTA reads (no separate sample tag FASTQ input); enable with `use_sampletags`/`skip_sampletags` and `species`.

## v0.2.0 - 2026-04-24

- Sample config vocabulary: replaced `bead_version` with `allowedlist` (96/384) and `diversity_insets` (yes/no). Both are optional; bead chemistry is auto-detected from R1 linkers and the declared fields are used as a QC check and for logging.
- Added per-sample `use_sampletags` (yes/no). Sampletag demultiplexing is gated per-sample; `species` is only required when `use_sampletags` is yes.
- Added `cb_umi_max_errors` config key (integer, default 0) exposing cutadapt `-e` for the R1 linker trim.
- Added optional paired-fastq `downsample` (percentage in (0, 100], default 100) with `downsample_seed`. Uses seqtk; per-sample `uses.downsample` overrides the global value.
- Added per-sample linker QC report (`{sample}_linker_qc.html` under `linker_qc/`). The existing bead-class scan is reused to emit a per-read hamming-distance histogram for both v1 and enhanced chemistries. Helps pick a value for `cb_umi_max_errors`.
- Added cross-aligner per-cell correlation distribution plot in the per-sample comparison report (Pearson on log1p counts, one density per aligner pair, matched by cell barcode).
- SBG SCE generator remaps rownames to Ensembl IDs against the STARsolo `features.tsv` and keeps the BD symbol in `rowData$name`. Marker-voting cell types are now assigned on SBG cells.
- HeLa use case (`use_case: hela`) added; cell-cycle phase panels replace marker cell-type panels on clonal samples.
- Paper figures: pairwise pseudobulk scatter matrices per dataset (`*_bio_pseudobulk_correlation.pdf`) replace the old Pearson-r heatmaps for HeLa and sendoel.
- Paper figures: cluster confusion matrices reordered by Hungarian matching (`clue::solve_LSAP`), 3x2 panel grid, per-tile counts, `log1p` fill.
- Per-cell cross-aligner Pearson r computed on sparse matrices; scales to the >40k-cell sendoel P60 sample.

## v0.1.0

First tagged release of rhapsodist, a Snakemake workflow for single-cell BD Rhapsody data. 

It runs starsolo, kallisto and alevin-fry on the same inputs (and optionally the BD Rhapsody pipeline in a Singularity container) from either SRA accessions or local fastq paths, standardises barcodes across aligners, and produces per-sample reports covering QC, pseudobulk agreement (MARD and correlation), per-barcode UMI concordance (Bland-Altman), barcode overlap, per-pipeline Seurat clustering and UMAPs, marker-based cell-type assignment, marker co-expression (UpSet and combination bars), marker dotplots, cluster and cell-type adjusted Rand indices, guide RNA join rates for pooled CRISPR screens, and wall-clock / peak-memory benchmarks per pipeline step. Includes a simulated dataset for CI and integration tests.

Rhapsodist can be used to benchmark/compare aligners/pipelines, or to run one of them.
