#!/bin/bash
## Simulated reads, with and without zero counts: sbatch slurm/02_simulations.sh
#SBATCH --job-name=rhapsodist-sim
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000
#SBATCH --tmp=20000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

run_workflow workflow/Snakefile configs/sim_config.yaml
run_workflow workflow/Snakefile configs/sim_sparse_config.yaml
