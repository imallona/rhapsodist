#!/bin/bash
## Both datasets at 1 and 2 linker mismatches, one array task per config.
## Uses the fastqs fetched by 03 and 04: sbatch slurm/05_linker_tolerance.sh
#SBATCH --job-name=rhapsodist-linker
#SBATCH --array=0-3
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=50000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%A_%a.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

configfiles=(
    configs/moro_mallona2025_linker1_config.yaml
    configs/moro_mallona2025_linker2_config.yaml
    configs/sendoel2024_linker1_config.yaml
    configs/sendoel2024_linker2_config.yaml
)

run_workflow workflow/Snakefile "${configfiles[$SLURM_ARRAY_TASK_ID]}"
