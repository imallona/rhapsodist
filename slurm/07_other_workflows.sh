#!/bin/bash
## UniverSC, zUMIs and OpenPipelines on 10% of the HeLa (task 0) and
## epidermis (task 1) reads, after 03 and 04: sbatch slurm/07_other_workflows.sh
## A tool without support for the bead version fails; the others still run.
#SBATCH --job-name=rhapsodist-others
#SBATCH --array=0-1
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8000
#SBATCH --tmp=200000
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%A_%a.out
#SBATCH --signal=B:TERM@300

source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/common.sh"

if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    configfile=configs/moro_mallona2025_config.yaml
    working_dir=output/moro_mallona2025
    sample=hela_unmod
    workflows=(universc openpipelines)
    settings=(
        universc_technology=bd-rhapsody-v2
        openpipelines_bead_version=EnhV2
        openpipelines_whitelist_dir=workflow/data/whitelist_384x3
    )
else
    configfile=configs/sendoel2024_downsampled_config.yaml
    working_dir=output/sendoel2024_downsampled
    sample=sample_16_wta_p60
    workflows=(universc zumis openpipelines)
    settings=(
        universc_technology=bd-rhapsody
        openpipelines_bead_version=Enh
    )
    ## The four aligners on the same 10% of the reads.
    run_workflow workflow/Snakefile "$configfile"
fi

targets=()
for workflow in "${workflows[@]}"; do
    targets+=("$working_dir/$workflow/$sample/${sample}_${workflow}_sce.rds")
done

run_workflow paper/Snakefile "$configfile" --keep-going "${targets[@]}" \
    --config universc_sandbox="$UNIVERSC_SANDBOX" zumis_dir="$ZUMIS_DIR" "${settings[@]}"
