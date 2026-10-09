#!/bin/bash
## Builds the conda environments and pulls the BD image, so the timed runs
## do not: sbatch slurm/01_environments.sh
#SBATCH --job-name=rhapsodist-envs
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000
#SBATCH --tmp=20000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

## A simulated, an SRA and a sample tag config.
for configfile in configs/sim_config.yaml configs/moro_mallona2025_config.yaml \
        configs/gse282765_config.yaml; do
    run_workflow workflow/Snakefile "$configfile" --conda-create-envs-only
done

## cwltool reads the image from this file in the repository root.
bd_image="bdgenomics_rhapsody:2.2.1.sif"
if [ ! -f "$bd_image" ]; then
    apptainer pull "$bd_image" docker://bdgenomics/rhapsody:2.2.1
fi
