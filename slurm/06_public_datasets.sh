#!/bin/bash
## GSE282765 and GSE301173, one array task each. Uses the mouse reference
## fetched by 04: sbatch slurm/06_public_datasets.sh
#SBATCH --job-name=rhapsodist-public
#SBATCH --array=0-1
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=200000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%A_%a.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

configfiles=(
    configs/gse282765_config.yaml
    configs/gse301173_config.yaml
)

run_workflow workflow/Snakefile "${configfiles[$SLURM_ARRAY_TASK_ID]}"
