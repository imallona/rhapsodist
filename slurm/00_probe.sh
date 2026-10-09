#!/bin/bash
## Checks the job environment before the runs: sbatch slurm/00_probe.sh
#SBATCH --job-name=rhapsodist-probe
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000
#SBATCH --tmp=5000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%j.out

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

echo "snakemake $(snakemake --version)"
echo "apptainer $(apptainer --version 2>&1)"
echo "nextflow $(nextflow -version 2>&1 | sed -n 's/.*version //p' | head -n 1)"
echo "TMPDIR ${TMPDIR:-unset}, $(df -h --output=avail "${TMPDIR:-/tmp}" | tail -n 1) free"

echo "== internet"
curl -fsSI https://ftp.ebi.ac.uk/pub/databases/gencode/ | head -n 1 \
    || echo "no internet access"

echo "== sandbox build"
apptainer build --sandbox "${TMPDIR:-/tmp}/probe_sandbox" docker://alpine:3.20 \
    && echo "sandbox build works" || echo "sandbox build failed"

echo "== dry run"
run_workflow workflow/Snakefile configs/sim_config.yaml --dry-run --quiet
