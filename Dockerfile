## rhapsodist with Snakemake and its pinned conda environments. The BD pipeline
## (sbg aligner) is not included: it starts its own containers.
FROM mambaorg/micromamba:2.9.0
ENV SNAKEMAKE_CONDA_PREFIX=/opt/rhapsodist_envs
USER root
COPY . /opt/rhapsodist
WORKDIR /opt/rhapsodist
RUN bash container_build.sh
WORKDIR /work
ENTRYPOINT ["micromamba", "run", "-n", "base", "rhapsodist"]
