#!/usr/bin/env bash
## Runs zUMIs on one sample of BD Rhapsody v1 reads (cell label at read 1
## positions 1-9, 22-30 and 44-52, UMI at 53-60) with the conda environment of
## ZUMIS_DIR, a clone: git clone --branch 2.9.7 https://github.com/sdparekh/zUMIs
##
## run_zumis.sh ZUMIS_DIR WORKDIR SAMPLE R1 R2 GENOME_FA GTF CORES MEM_GB
##
## Builds a STAR index with zUMIs's STAR in WORKDIR on first use. Removes the
## fastq header comments first: zUMIs keeps a comment that contains a space in
## the read name, and STAR then fails. Counts: WORKDIR/zUMIs_output/expression.
set -euo pipefail

if [ $# -ne 9 ]; then
    sed -n '6p' "$0" >&2
    exit 1
fi
zumis_dir=$(realpath "$1"); workdir=$(realpath -m "$2"); sample=$3
r1=$(realpath "$4"); r2=$(realpath "$5"); genome=$(realpath "$6"); gtf=$(realpath "$7")
cores=$8; mem_gb=$9

if [ ! -f "$zumis_dir/zUMIs.sh" ]; then
    echo "$zumis_dir is not a zUMIs clone" >&2
    exit 1
fi
mkdir -p "$workdir"
cd "$workdir"

if [ ! -d "$zumis_dir/zUMIs-env/bin" ]; then
    mkdir -p "$zumis_dir/zUMIs-env"
    cat "$zumis_dir"/zUMIs-miniconda.parta* | tar -xj -C "$zumis_dir/zUMIs-env"
fi
export PATH="$zumis_dir/zUMIs-env/bin:$PATH"

cdna_length=$(zcat "$r2" 2> /dev/null | awk 'NR == 2 {print length($0); exit}' || true)

if [ ! -f star_index/SA ]; then
    mkdir -p star_index
    STAR --runMode genomeGenerate --runThreadN "$cores" --genomeDir star_index \
        --genomeFastaFiles "$genome" --sjdbGTFfile "$gtf" \
        --sjdbOverhang $((cdna_length - 1)) \
        --limitGenomeGenerateRAM $((mem_gb * 1000000000)) > star_index.log 2>&1
fi

strip_header_comments() {
    zcat "$1" | awk 'NR % 4 == 1 {print $1; next} {print}' | pigz -p "$cores" > "$2"
}
strip_header_comments "$r1" "${sample}_R1.fastq.gz"
strip_header_comments "$r2" "${sample}_R2.fastq.gz"

cat > "${sample}.yaml" <<YAML
project: $sample
sequence_files:
  file1:
    name: $workdir/${sample}_R1.fastq.gz
    base_definition:
      - BC(1-9,22-30,44-52)
      - UMI(53-60)
  file2:
    name: $workdir/${sample}_R2.fastq.gz
    base_definition:
      - cDNA(1-$cdna_length)
reference:
  STAR_index: $workdir/star_index
  GTF_file: $gtf
  exon_extension: no
  extension_length: 0
  scaffold_length_min: 0
  additional_files:
  additional_STAR_params:
out_dir: $workdir
num_threads: $cores
mem_limit: $mem_gb
filter_cutoffs:
  BC_filter:
    num_bases: 1
    phred: 20
  UMI_filter:
    num_bases: 1
    phred: 20
barcodes:
  barcode_num: null
  barcode_file: null
  barcode_sharing: null
  automatic: yes
  BarcodeBinning: 1
  nReadsperCell: 100
  demultiplex: no
counting_opts:
  introns: yes
  intronProb: no
  downsampling: 0
  strand: 1
  Ham_Dist: 0
  velocyto: no
  primaryHit: yes
  multi_overlap: no
  fraction_overlap: 0
  twoPass: yes
make_stats: yes
which_Stage: Filtering
YAML

bash "$zumis_dir/zUMIs.sh" -c -y "${sample}.yaml" > "${sample}_zumis.log" 2>&1 || true
if [ ! -s "zUMIs_output/expression/${sample}.dgecounts.rds" ]; then
    echo "zUMIs wrote no count table, see $workdir/${sample}_zumis.log" >&2
    exit 1
fi
