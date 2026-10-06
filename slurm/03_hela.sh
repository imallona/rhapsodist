#!/bin/bash
## HeLa, 10% of the reads, four aligners: sbatch slurm/03_hela.sh
#SBATCH --job-name=rhapsodist-hela
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=200000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

run_workflow workflow/Snakefile configs/moro_mallona2025_config.yaml
