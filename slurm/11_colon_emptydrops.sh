#!/bin/bash
## GSE282765 with emptyDrops for every aligner, in its own working directory,
## so the aligners are compared on cells called the same way. Writes under
## RUN_ROOT when site.env sets it, else under output/, and on success copies
## the results to RESULTS_DIR:
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

working_dir="${RUN_ROOT:-output}/gse282765_emptydrops"
run_workflow workflow/Snakefile configs/gse282765_config.yaml \
    --config cell_filtering=emptydrops working_dir="$working_dir"
archive_results "$working_dir" "$RESULTS_DIR/$(git rev-parse --short HEAD)/gse282765_emptydrops"
