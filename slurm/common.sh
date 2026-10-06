## Shared setup of the job scripts. Sourced from the repository root.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

if [ -f slurm/site.env ]; then
    source slurm/site.env
fi

if [ -n "${PROXY_MODULE:-}" ]; then
    module load "$PROXY_MODULE"
fi

if ! command -v snakemake > /dev/null && [ -n "${CONDA_INIT:-}" ]; then
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
fi

CORES="${SLURM_CPUS_PER_TASK:-10}"

conda_prefix_flag=()
if [ -n "${CONDA_PREFIX_DIR:-}" ]; then
    mkdir -p "$CONDA_PREFIX_DIR"
    conda_prefix_flag=(--conda-prefix "$CONDA_PREFIX_DIR")
fi

if [ -n "${APPTAINER_CACHEDIR:-}" ]; then
    mkdir -p "$APPTAINER_CACHEDIR"
    export APPTAINER_CACHEDIR
    export NXF_SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
fi

echo "host $(hostname), job ${SLURM_JOB_ID:-none}, commit $(git rev-parse --short HEAD)"
echo "cpu $(sed -n 's/^model name[^:]*: //p' /proc/cpuinfo | head -n 1), cores $CORES"

## run_workflow SNAKEFILE CONFIGFILE [snakemake arguments]
## Slurm sends TERM before the time limit; it is forwarded to Snakemake.
run_workflow() {
    local snakefile=$1 configfile=$2
    shift 2
    echo "$(date -Is) $snakefile $configfile $*"
    snakemake --snakefile "$snakefile" --configfile "$configfile" \
        --use-conda "${conda_prefix_flag[@]}" --benchmark-extended \
        --cores "$CORES" --rerun-incomplete "$@" &
    local child=$! status=0 stopping=""
    trap 'stopping=yes; kill -TERM "$child" 2> /dev/null' TERM
    wait "$child" || status=$?
    if [ -n "$stopping" ]; then
        wait "$child" || status=$?
    fi
    trap - TERM
    return "$status"
}
