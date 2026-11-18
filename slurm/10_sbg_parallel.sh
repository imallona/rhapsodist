#!/bin/bash
## Reruns only the BD pipeline of each dataset, to time it with
## cwl-runner --parallel. Fetches reads and builds the reference again when the
## scratch purge removed them: sbatch slurm/10_sbg_parallel.sh
#SBATCH --job-name=rhapsodist-sbg-parallel
#SBATCH --array=0-6
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=200000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%A_%a.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

configfiles=(
    configs/sim_config.yaml
    configs/sim_sparse_config.yaml
    configs/moro_mallona2025_config.yaml
    configs/sendoel2024_config.yaml
    configs/sendoel2024_downsampled_config.yaml
    configs/gse282765_config.yaml
    configs/gse301173_config.yaml
)
samples=(simulated simulated hela_unmod sample_16_wta_p60 sample_16_wta_p60 colon neutrophils)

rule="run_sbg_cwl_${samples[$SLURM_ARRAY_TASK_ID]}"
run_workflow workflow/Snakefile "${configfiles[$SLURM_ARRAY_TASK_ID]}" \
    --forcerun "$rule" --until "$rule"
