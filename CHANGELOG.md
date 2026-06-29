# Changelog

## v0.3.0 - 2026-06-29

### Requests

- Added `alevin_usa`: alevin-fry quantifies in USA mode and adds `spliced`, `unspliced` and `ambiguous` assays. Requires `alevin_sketch`.
- Documented that BD sample tags come from the WTA reads, with no separate tag FASTQ input.
- The transcriptome input now accepts a plain `.fa`, not only `.fa.gz`.
- A sample can list several fastqs per mate; they are concatenated before processing, paired by order.
- Sample tags can be called without starsolo, via an alignment-free search that writes the same count table.
- The cross-pipeline comparison report is built only with two or more aligners; single-aligner runs skip it.
- Fixed sample tag alignment that rejected every read as too short; STAR now filters on matched bases.
- Added per-sample `sampletags` to declare which tags a sample carries, as a list or tag-to-label mapping.
- Each aligner's counts are split into one HDF5 SCE per tag; demultiplexing moved from the report into `demux_sampletags.R`.
- Downgraded the setuptools pin in the pyroe env to fix a `pkg_resources` `ModuleNotFoundError`.
- Sped up the per-sampletag split with `sampletag_split_backend: memory | delayed` (default `memory`); memory reads the source once.
- Added `output_format: sce | h5ad | both` (default `sce`); h5ad written via anndataR next to each SCE and split.

### Extras

- Extended the integration tests: memory and delayed sampletag scenarios, dry-run coverage, an `integration` PR label, and failure artifact uploads.
- Fixed the simulated sampletag path: `sampletag_fa` resolved from the wrong directory, so no sampletag reads were produced.
- Pinned the self-compiled kallisto to `v0.52.0` and bustools to `0.45.1` from the kallisto conda env.
- `get_txp2gene` reads `transcript_id` and `gene_id` GTF attributes by name, so any attribute order works.
- Added `cell_filtering: none` to keep every observed barcode across STARsolo, alevin and kallisto.

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
