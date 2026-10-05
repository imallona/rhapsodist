#!/usr/bin/env bash
## Runs UniverSC on one sample from a writable sandbox of its image, which
## UniverSC needs to rewrite the Cell Ranger whitelist:
## apptainer build --sandbox SANDBOX docker://tomkellygenetics/universc:1.2.7
##
## run_universc.sh SANDBOX WORKDIR SAMPLE TECHNOLOGY R1 R2 GENOME_FA GTF CORES MEM_GB
##
## Builds the Cell Ranger reference in WORKDIR on first use. If Cell Ranger
## fails, pairs whose converted read 1 has sequence and quality of unequal
## length are dropped and Cell Ranger runs again; their number goes to
## WORKDIR/SAMPLE_dropped_pairs.txt. Counts: WORKDIR/SAMPLE_outs.
set -euo pipefail

if [ $# -ne 10 ]; then
    sed -n '6p' "$0" >&2
    exit 1
fi
sandbox=$1; workdir=$(realpath "$2"); sample=$3; technology=$4
r1=$(realpath "$5"); r2=$(realpath "$6"); genome=$(realpath "$7"); gtf=$(realpath "$8")
cores=$9; mem_gb=${10}

if [ ! -d "$sandbox/universc" ]; then
    echo "$sandbox is not a UniverSC sandbox" >&2
    exit 1
fi

mkdir -p "$workdir/${sample}_fastq" "$sandbox/work" "$sandbox/reads" "$sandbox/reference"
ln -sf "/reads/$(basename "$r1")" "$workdir/${sample}_fastq/${sample}_S1_L001_R1_001.fastq.gz"
ln -sf "/reads/$(basename "$r2")" "$workdir/${sample}_fastq/${sample}_S1_L001_R2_001.fastq.gz"

apptainer exec --writable \
    -B "$workdir:/work" -B "$(dirname "$r1"):/reads" -B "$(dirname "$genome"):/reference" \
    "$sandbox" bash -s "$sample" "$technology" "$(basename "$genome")" "$(basename "$gtf")" \
    "$cores" "$mem_gb" <<'IN_CONTAINER'
set -euo pipefail
sample=$1; technology=$2; genome=$3; gtf=$4; cores=$5; mem_gb=$6
cd /work

if [ ! -f reference/reference.json ]; then
    cellranger mkref --genome=reference --fasta="/reference/$genome" \
        --genes="/reference/$gtf" --nthreads="$cores" --memgb="$mem_gb" > mkref.log 2>&1
fi

bash /universc/launch_universc.sh -t "$technology" \
    -R1 "${sample}_fastq/${sample}_S1_L001_R1_001.fastq.gz" \
    -R2 "${sample}_fastq/${sample}_S1_L001_R2_001.fastq.gz" \
    -i "$sample" -r reference --jobmode local \
    --localcores "$cores" --localmem "$mem_gb" > "${sample}_universc.log" 2>&1 || true

if [ -f "$sample/outs/filtered_feature_bc_matrix/matrix.mtx.gz" ]; then
    ln -sfn "$sample/outs" "${sample}_outs"
    exit 0
fi

converted="input4cellranger_${sample}"
fixed="${sample}_fixed_fastq"
mkdir -p "$fixed"
paste <(paste - - - - < "$converted/${sample}_S1_L001_R1_001.fastq") \
      <(paste - - - - < "$converted/${sample}_S1_L001_R2_001.fastq") |
    awk -F'\t' -v r1="$fixed/${sample}_S1_L001_R1_001.fastq" \
               -v r2="$fixed/${sample}_S1_L001_R2_001.fastq" \
               -v dropped_file="${sample}_dropped_pairs.txt" '
        length($2) == length($4) {
            print $1 "\n" $2 "\n" $3 "\n" $4 > r1
            print $5 "\n" $6 "\n" $7 "\n" $8 > r2
            next
        }
        { dropped++ }
        END { print dropped + 0 > dropped_file }'
gzip -f -1 "$fixed"/*.fastq

cellranger count --id="${sample}_fixed" --fastqs="$fixed" --sample="$sample" \
    --transcriptome=reference --chemistry=SC3Pv2 \
    --localcores="$cores" --localmem="$mem_gb" > "${sample}_cellranger_fixed.log" 2>&1
ln -sfn "${sample}_fixed/outs" "${sample}_outs"
IN_CONTAINER
