#!/bin/bash
## Reruns only the BD pipeline of each dataset, to time it with
## cwl-runner --parallel. Fetches reads and builds the reference again when
## they are missing: sbatch slurm/10_sbg_parallel.sh
## With RUN_ROOT set in site.env the run writes there, and the BD timings are
## copied back to the benchmarks of the config's working_dir.
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

configfile=${configfiles[$SLURM_ARRAY_TASK_ID]}
sample=${samples[$SLURM_ARRAY_TASK_ID]}
rule="run_sbg_cwl_$sample"
working_dir=$(working_dir_of "$configfile")
run_workflow workflow/Snakefile "$configfile" \
    --config working_dir="$working_dir" --forcerun "$rule" --until "$rule"

if [ -n "${RUN_ROOT:-}" ]; then
    benchmarks="$(sed -n 's/^working_dir: *//p' "$configfile")/benchmarks"
    mkdir -p "$benchmarks"
    cp -v "$working_dir"/benchmarks/sbg_cwl_{,steps_}"$sample".* "$benchmarks"/
fi
