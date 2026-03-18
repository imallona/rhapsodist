# Aim

Rhapsodist is a Snakemake workflow to process BD Rhapsody WTA (enhanced beads) data.

## TL/DR

Simulations:

```
snakemake --use-conda --cores 10 --configfile sim_config.yaml
```

Real data:

```
snakemake --use-conda --cores 10 --configfile config.yaml
```

## Configuration

Copy `config.yaml` and fill in the fields below before running on real data.

**Paths and resources**

- `Rbin`: path to the R binary (e.g. `/usr/bin/R`)
- `STAR`: path to the STAR binary, or `STAR` if on `PATH`
- `nthreads`: number of CPU threads
- `max_mem_mb`: RAM limit in MB
- `working_dir`: absolute path where outputs will be written

**Reference files** (uncompressed unless noted)

- `gtf_origin`: `"gencode"` or `"ensembl"`
- `gtf`: path to GTF annotation file (uncompressed)
- `genome`: path to genome FASTA (uncompressed)
- `transcriptome`: path to transcriptome FASTA (can be gzipped)
- `sjdbOverhang`: read length minus 1 (e.g. 70 for 71 bp reads)

**Aligners**

- `aligner`: list of one or more aligners to run, e.g. `['starsolo', 'kallisto', 'alevin']`

**Samples**

```yaml
samples:
  - name: my_sample
    uses:
      cb_umi_fq: /path/to/R1.fastq.gz   # barcode + UMI read
      cdna_fq: /path/to/R2.fastq.gz     # cDNA read
      whitelist: 384x3                   # 384x3 for Enhanced beads, 96x3 for v1/Enh beads
      species: human                     # human or mouse
```

**BD Rhapsody official pipeline (optional)**

To also run the SBG/CWL pipeline, set `sbg_cwl` to the CWL workflow file path. For real data, provide either `sbg_reference_url` (URL to download the reference archive) or `sbg_reference_archive` (path to a local copy). For ingest mode (skip re-running CWL), set `sbg_mex_dir` to a directory of pre-computed MEX output.

## Contributors

- Izaskun Mallona
- Jiayi Wang
- Giulia Moro

Tools used include STAR, subread (featureCounts), samtools, kallisto, and alevin.

## Contact

izaskun.mallona at mls.uzh.ch, Mark D. Robinson lab
https://www.mls.uzh.ch/en/research/robinson.html

## Started

30 July 2024, keeping history from https://github.com/imallona/rock_roi_method
