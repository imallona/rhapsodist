# paper/

Manuscript-only scaffolding. Reads RDS outputs produced by the main rhapsodist pipeline and writes CSV and PDF artefacts ready for inclusion in the paper. Nothing in `workflow/` depends on this directory.

## Layout

- `Snakefile`: optional, three rules (`simulations_figure`, `sendoel_figure`, `benchmarks_figure`).
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
