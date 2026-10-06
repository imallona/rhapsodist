#!/bin/bash
## Figure tables and panels of the simulated, HeLa and epidermis runs:
## sbatch slurm/08_figures.sh
#SBATCH --job-name=rhapsodist-figures
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

for configfile in configs/sim_config.yaml configs/sim_sparse_config.yaml \
        configs/moro_mallona2025_config.yaml configs/sendoel2024_config.yaml; do
    run_workflow paper/Snakefile "$configfile"
done
