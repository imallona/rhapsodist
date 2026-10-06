#!/bin/bash
## Copies reports, count objects, benchmarks, figures and logs to
## RESULTS_DIR/<commit>: sbatch slurm/09_archive.sh
#SBATCH --job-name=rhapsodist-archive
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

destination="$RESULTS_DIR/$(git rev-parse --short HEAD)"
mkdir -p "$destination"

rsync -a --prune-empty-dirs \
    --include='*/' --include='*.html' --include='*_sce.rds' --include='*.h5ad' \
    --include='run_info.tsv' --include='benchmarks/***' --include='paper/***' \
    --include='logs/***' --exclude='*' \
    output/ "$destination/"
rsync -a slurm/logs/ "$destination/slurm_logs/"
