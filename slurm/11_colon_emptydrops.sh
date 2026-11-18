#!/bin/bash
## GSE282765 with emptyDrops for every aligner, in its own working directory,
## so the aligners are compared on cells called the same way:
## sbatch slurm/11_colon_emptydrops.sh
#SBATCH --job-name=rhapsodist-colon-emptydrops
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=200000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

run_workflow workflow/Snakefile configs/gse282765_config.yaml \
    --config cell_filtering=emptydrops working_dir=output/gse282765_emptydrops
