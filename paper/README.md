# paper/

Manuscript-only scaffolding. Reads RDS outputs produced by the main rhapsodist pipeline and writes CSV and PDF artefacts ready for inclusion in the paper. Nothing in `workflow/` depends on this directory.

## Layout

- `Snakefile`: the figure rules (`simulations_figure`, `experimental_figure`, `linker_qc_figure`, `benchmarks_figure`), the comparison rules (`compare_linker_tolerance`, `compare_workflows`, `compare_sampletag_calls`, each on when the config sets `linker_tolerance_runs`, `other_workflows` or `deposited_sampletag_files`), and the other workflow rules.
- `assemble_paper_figures.R`: the single script that does the work.
- `export_helpers.R`: small helpers for CSV and PDF writing, bootstrap confidence intervals, and ggrastr-based rasterisation.
- `env.yaml`: conda environment used by the Snakemake rules. Lighter than `workflow/envs/r_bioc.yaml` and does not need to track it.

## Usage

Run it after the main pipeline has completed, passing the same config. All commands below must be run from the rhapsodist repo root, not from inside `paper/`, because the configs use relative paths such as `working_dir: output/simul` that snakemake resolves against the current working directory.

```bash
cd /home/imallona/src/rhapsodist   # repo root, not paper/
source ~/miniconda3/bin/activate
conda activate snakemake

# Simulation dataset: working_dir is output/simul, sample is "simulated".
# Outputs land under output/simul/paper/simulated/ and output/simul/paper/.
snakemake -s paper/Snakefile \
    --configfile configs/sim_config.yaml \
    --use-conda --cores 4

# Real dataset: working_dir is output/sendoel2024, sample is
# "sample_16_wta_p60". Outputs land under
# output/sendoel2024/paper/sample_16_wta_p60/ and
# output/sendoel2024/paper/.
snakemake -s paper/Snakefile \
    --configfile configs/sendoel2024_config.yaml \
    --use-conda --cores 4
```

In general, per-sample materials land under `<working_dir>/paper/<sample>/` and shared benchmark tables under `<working_dir>/paper/`, where `<working_dir>` is whatever the config sets (relative paths resolve from the repo root).

## Running the script directly

The Snakefile is just a thin wrapper. You can call the script manually, again from the repo root:

```bash
cd /home/imallona/src/rhapsodist
Rscript paper/assemble_paper_figures.R \
    --kind simulation \
    --working_dir output/simul \
    --sample simulated \
    --aligners starsolo,kallisto,alevin,sbg \
    --n_expected_cells 1000
```

Three kinds are supported: `simulation`, `biology`, `benchmarks`.

## Inputs it expects

The script reads only files the main pipeline writes:

- SCE RDS files at `<wd>/<aligner>/<sample>/<sample>_<aligner>_sce.rds`
- Simulation truth at `<wd>/simulate/cell_barcodes.txt` and `<wd>/simulate/true_mex/` (matrix market format)
- Biology derived objects at `<wd>/<sample>_biology_*.rds` including `qc`, `cb_umi`, `clusters`, `pseudobulk_mard`, and per-pipeline Seurat caches `<sample>_biology_<aligner>_seurat.rds`
- Benchmark files under `<wd>/benchmarks/*.txt`
- The machine record at `<wd>/run_info.tsv`; the `benchmarks` kind copies it to `paper/run_info.csv` and writes `bench_stage.csv` from `workflow/data/benchmark_stages.tsv`; report and install rules are left out of the time and memory tables

## Sample tag concordance

`compare_sampletag_calls.R` compares the sample tag calls of a run with per-tag expression files of the BD pipeline, such as the ones deposited with GSM8696920. It writes a CSV and a PDF with, per deposited tag, the cells called with the same tag, the other tag, multiplet, undetermined or absent.

```bash
Rscript paper/compare_sampletag_calls.R \
    --demux output/gse282765/sampletags/colon/sampletag_demux.tsv.gz \
    --deposited mouse_sampletag_3=GSM8696920_Colon_healthy_Mm_PHIL_ST03_Expression_Data.st.txt.gz \
                mouse_sampletag_4=GSM8696920_Colon_healthy_Mm_WT_ST04_Expression_Data.st.txt.gz \
    --whitelist_dir workflow/data/whitelist_384x3 \
    --out_prefix output/gse282765/paper/colon/sampletag_concordance
```

## Linker tolerance

`compare_linker_tolerance.R` compares runs of one sample that differ in `cb_umi_max_errors`. Per aligner and setting it writes cells, median UMIs per cell, and Spearman and scaled MARD of the pseudobulk against the first run, as a CSV and a PDF. With `--reference_aligner sbg` it also writes `<out_prefix>_sbg.csv`: Spearman, scaled MARD and MARD of every aligner and setting against `sbg` of the first run, on gene symbols shared by all of them. The paper rule passes it when the run has `sbg`.

```bash
Rscript paper/compare_linker_tolerance.R \
    --runs 0=output/sendoel2024 1=output/sendoel2024_linker1 2=output/sendoel2024_linker2 \
    --sample sample_16_wta_p60 --aligners starsolo,kallisto,alevin \
    --reference_aligner sbg \
    --out_prefix output/sendoel2024/paper/sample_16_wta_p60/linker_tolerance
```

## Workflow comparison

`compare_workflows.R` compares UniverSC, zUMIs and OpenPipelines with the rhapsodist aligners on one sample. Per workflow and aligner it writes cells, cells shared by barcode, Spearman and scaled MARD of the pseudobulk, and the time and memory of the workflow's benchmark file. A workflow without an SCE is listed as failed.

```bash
Rscript paper/compare_workflows.R \
    --working_dir output/sendoel2024_downsampled40 --sample sample_16_wta_p60 --bead_version v1 \
    --aligners starsolo,kallisto,alevin,sbg --workflows universc,zumis,openpipelines \
    --out_prefix output/sendoel2024_downsampled40/paper/sample_16_wta_p60/workflow_comparison
```

## UniverSC

UniverSC rewrites the Cell Ranger whitelist inside its image, so `universc/run_universc.sh` runs it from a writable copy:

```bash
apptainer build --sandbox universc_sandbox docker://tomkellygenetics/universc:1.2.7
apptainer exec --writable universc_sandbox python2.7 -m pip install scipy==1.2.3
```

- scipy: the image has 0.18.1, and the `FILTER_BARCODES` stage of Cell Ranger fails with it on the epidermis sample. The sandbox gets 1.2.3, the last release for Python 2.7.
- Config: `universc_sandbox` and `universc_technology` (`bd-rhapsody` or `bd-rhapsody-v2`).
- Target: `<working_dir>/universc/<sample>/<sample>_universc_sce.rds`, from the downsampled reads.
- Enhanced beads: UniverSC 1.2.7 cannot run them.
  - `launch_universc.sh` line 3175 tests for `bd-rhapsody`, so `-t bd-rhapsody-v2` removes no linker.
  - Its pattern needs two bases before the first cell label section: 22% of HeLa read 1 records match, those with the `GT` inset.
  - `whitelists/bd_rhapsody_cell_label_section*.txt` hold 96 sequences per section, not 384.
  - With line 3175 corrected and the 384 lists, Cell Ranger finds 0.00% valid barcodes on 2 M HeLa read pairs.
- v1 beads: the conversion writes some read 1 records with sequence and quality of unequal length, and Cell Ranger fails. The script drops those pairs, counts them in `<sample>_dropped_pairs.txt` and reruns Cell Ranger.

## zUMIs

`zumis/run_zumis.sh` runs zUMIs 2.9.7 from a clone, with the conda environment it ships:

```bash
git clone --branch 2.9.7 https://github.com/sdparekh/zUMIs
```

- Config: `zumis_dir`, the clone.
- Target: `<working_dir>/zumis/<sample>/<sample>_zumis_sce.rds`, from the downsampled reads.
- v1 beads only: zUMIs reads the cell label at fixed positions, which the variable-length inset of Enhanced beads shifts.
- Fastq header comments are removed first. zUMIs keeps a comment that contains a space, as in SRA headers, in the read name, and STAR then fails.

## OpenPipelines

`openpipelines/run_openpipelines.sh` runs the `bd_rhapsody` component of OpenPipelines v4.2.0 with nextflow and apptainer. The component runs the BD Rhapsody Sequence Analysis CWL pipeline v2.2.1 in the image `bdgenomics/rhapsody:2.2.1`, the pipeline of the `sbg` aligner.

- Config: `openpipelines_reference` (default: the `sbg` reference archive), `openpipelines_bead_version` (`Enh` or `EnhV2`), `openpipelines_whitelist_dir` for `EnhV2`.
- Target: `<working_dir>/openpipelines/<sample>/<sample>_openpipelines_sce.rds`, from the reads of the `sbg` aligner.
- On the simulated reads its counts equal those of the `sbg` aligner.

## Manuscript figures

`manuscript_figures.R` draws the figures of the paper at print size (6.3 in wide, 7 pt text) from CSV files only. It expects one folder per run under `--results`, named as the working directories (`simul`, `simul_sparse`, `moro_mallona2025`, `sendoel2024`, `sendoel2024_downsampled40`, `gse282765`, `gse301173`), each with its `paper/` folder. Two scripts write the CSVs that the other rules do not:

- `export_embeddings.R` writes `bio_umap.csv`: UMAP coordinates, cluster, cell type and cell cycle phase per cell and aligner, from the Seurat objects and clusters table of `04_biology.Rmd`. It needs Seurat, so run it in the `workflow/envs/r_bioc.yaml` environment.
- `export_pseudobulk.R` writes `bio_pseudobulk.csv`, the summed counts per aligner on shared gene symbols, and with `--cells_out` the cells called per aligner.

```bash
Rscript paper/export_embeddings.R --working_dir output/sendoel2024 \
    --sample sample_16_wta_p60 --out output/sendoel2024/paper/sample_16_wta_p60/bio_umap.csv
Rscript paper/export_pseudobulk.R --working_dir output/gse282765 \
    --sample colon --out output/gse282765/paper/colon/bio_pseudobulk.csv
Rscript paper/manuscript_figures.R --results output --out figures
```

`--only fig2,linker_tolerance` draws a subset; the names are listed at the end of the script.
