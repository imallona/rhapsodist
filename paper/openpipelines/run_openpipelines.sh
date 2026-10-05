#!/usr/bin/env bash
## Runs the OpenPipelines v4.2.0 bd_rhapsody component on one sample, with
## nextflow and apptainer. The component runs the BD Rhapsody Sequence Analysis
## CWL pipeline v2.2.1 in the image bdgenomics/rhapsody:2.2.1.
##
## run_openpipelines.sh WORKDIR SAMPLE R1 R2 REFERENCE_ARCHIVE CORES MEM_GB
##
## The image is kept in NXF_SINGULARITY_CACHEDIR, or WORKDIR/images when unset.
## Outputs: WORKDIR/output/SAMPLE, the BD output directory, and WORKDIR/MEX,
## its filtered matrix unzipped, or the unfiltered one when BD calls no cells.
set -euo pipefail

if [ $# -ne 7 ]; then
    sed -n '6p' "$0" >&2
    exit 1
fi
workdir=$(realpath -m "$1"); sample=$2
r1=$(realpath "$3"); r2=$(realpath "$4"); reference=$(realpath "$5")
cores=$6; mem_gb=$7

mkdir -p "$workdir"
cd "$workdir"

cat > limits.config <<CONFIG
process.resourceLimits = [cpus: $cores, memory: ${mem_gb}.GB]
singularity.cacheDir = "${NXF_SINGULARITY_CACHEDIR:-$workdir/images}"
CONFIG

nextflow run openpipelines-bio/openpipeline -r v4.2.0 \
    -main-script target/nextflow/mapping/bd_rhapsody/main.nf \
    -profile singularity -c limits.config -work-dir "$workdir/nextflow_work" -resume \
    --id "$sample" --reads "$r1;$r2" --reference_archive "$reference" \
    --long_reads false --output_dir "$sample" --publish_dir "$workdir/output"

mex_zip=$(find "$workdir/output/$sample" -name '*_RSEC_MolsPerCell_MEX.zip' | head -n 1)
if [ -z "$mex_zip" ]; then
    mex_zip=$(find "$workdir/output/$sample" -name '*_RSEC_MolsPerCell_Unfiltered_MEX.zip' | head -n 1)
fi
basename "$mex_zip" > "$workdir/mex_source.txt"
mkdir -p "$workdir/MEX"
unzip -o "$mex_zip" -d "$workdir/MEX"
