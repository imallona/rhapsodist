#!/usr/bin/env bash
## Installs Snakemake and rhapsodist into the base environment of a micromamba
## image and builds the pinned conda environments of workflow/envs. Shared by
## Dockerfile and rhapsodist.def; run from the repository copy in the image.
set -euo pipefail

micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.12 snakemake=9.27.0 conda pip git
micromamba run -n base pip install --no-deps -e .

## the two configs together use every environment except sbg
for config in configs/sim_config.yaml configs/sendoel2024_config.yaml; do
    micromamba run -n base snakemake --snakefile workflow/Snakefile \
        --configfile "$config" --config "aligner=['starsolo','kallisto','alevin']" \
        --use-conda --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
        --conda-create-envs-only --cores 1
done
micromamba clean -a -y
