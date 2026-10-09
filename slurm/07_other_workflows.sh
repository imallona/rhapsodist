#!/bin/bash
## OpenPipelines on 10% of the HeLa reads (task 0), and UniverSC, zUMIs and
## OpenPipelines on 40% of the epidermis reads (task 1), after 03 and 04:
## sbatch slurm/07_other_workflows.sh
## UniverSC 1.2.7 has no whitelist for the 384x3 beads of the HeLa run; the
## comparison lists it as failed.
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
    workflows=(openpipelines)
    settings=(
        openpipelines_bead_version=EnhV2
        openpipelines_whitelist_dir=workflow/data/whitelist_384x3
    )
else
    configfile=configs/sendoel2024_downsampled_config.yaml
    working_dir=output/sendoel2024_downsampled40
    sample=sample_16_wta_p60
    workflows=(universc zumis openpipelines)
    settings=(
        universc_technology=bd-rhapsody
        openpipelines_bead_version=Enh
    )
    ## The four aligners on the same 40% of the reads.
    run_workflow workflow/Snakefile "$configfile"
fi

targets=()
for workflow in "${workflows[@]}"; do
    targets+=("$working_dir/$workflow/$sample/${sample}_${workflow}_sce.rds")
done

run_workflow paper/Snakefile "$configfile" --keep-going "${targets[@]}" \
    --config universc_sandbox="$UNIVERSC_SANDBOX" zumis_dir="$ZUMIS_DIR" "${settings[@]}"
