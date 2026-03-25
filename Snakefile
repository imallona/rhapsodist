#!/usr/bin/env snakemake -s
##
## Snakefile to process BD Rhapsody WTA data
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

import os.path as op
import os

configfile: "config.yaml"

## whitelists symlinking requires an absolute path before any include uses it
if not op.isabs(config['repo_path']):
    config['repo_path'] = op.join(workflow.basedir, config['repo_path'])

if not op.isabs(config['working_dir']):
    config['working_dir'] = op.join(workflow.basedir, config['working_dir'])

## when use_simulated is true, point genome/gtf/transcriptome/fastqs at generated outputs
if config.get('use_simulated', False):
    _sim_dir = op.join(config['working_dir'], 'simulate')
    config['genome']        = op.join(_sim_dir, 'genome.fa')
    config['gtf']           = op.join(_sim_dir, 'genes.gtf')
    config['transcriptome'] = op.join(_sim_dir, 'transcriptome.fa.gz')
    for _s in config['samples']:
        if not _s['uses'].get('cb_umi_fq'):
            _s['uses']['cb_umi_fq'] = op.join(_sim_dir, 'sim_R1.fq.gz')
        if not _s['uses'].get('cdna_fq'):
            _s['uses']['cdna_fq'] = op.join(_sim_dir, 'sim_R2.fq.gz')

include: "src/workflow_functions.py"

include: op.join('src', 'simulate.snmk')

try:
    os.makedirs(op.join(config['working_dir'], 'logs'), exist_ok=True)
    os.makedirs(op.join(config['working_dir'], 'benchmarks'), exist_ok=True)
except OSError as e:
    raise RuntimeError(
        f"Failed to create required directories under working_dir "
        f"'{config['working_dir']}': {e}"
    ) from e

print(get_sample_names())

_extra_targets = (
    [op.join(config['working_dir'], 'simulation_validation.html')]
    if config.get('use_simulated', False) else []
)

_sampletag_reports = (
    [] if config.get('skip_sampletags', False) else
    expand(op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_report.html'),
           sample = get_sample_names())
)

_has_sbg = 'sbg' in get_aligners()
_sbg_samples = get_sample_names() if _has_sbg else []

if _has_sbg and not config.get('sbg_cwl'):
    raise ValueError("'sbg' is in the aligner list but 'sbg_cwl' is not configured.")

_sbg_sce_targets = (
    expand(op.join(config['working_dir'], 'sbg', '{sample}', '{sample}_sbg_sce.rds'),
           sample = _sbg_samples)
    if _has_sbg else []
)

_comparison_reports = expand(
    op.join(config['working_dir'], '{sample}_comparison.html'),
    sample = get_sample_names()
)

_benchmarks_report = [op.join(config['working_dir'], 'benchmarks_report.html')]

rule all:
    input:
        expand(op.join(config['working_dir'], '{aligner}',  '{sample}', 'descriptive_report.html'),
               aligner = get_aligners(),
               sample = get_sample_names()),
        _sampletag_reports,
        _extra_targets,
        _sbg_sce_targets,
        _comparison_reports,
        _benchmarks_report
        # op.join(config['working_dir'], 'data', 'index', 'salmon', 'seq.bin'),
        # expand(op.join(config['working_dir'], 'alevin', '{sample}', 'alevin', 'quants_mat.gz'),
        #        sample = get_sample_names()),
        # expand(op.join(config['working_dir'], 'bustools', '{sample}', 'output.mtx'),
               # sample = get_sample_names()),
        # expand(op.join(config['working_dir'], 'rustody', '{sample}', 'flag'),
        #        sample = get_sample_names()),
        # expand(op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_counts.tsv.gz'),
        #        sample = get_sample_names()),
        # expand(op.join(config['working_dir'], 'starsolo', '{sample}', '{sample}_kallisto_sce.rds'),
        #        sample = get_sample_names()),
        # expand(op.join(config['working_dir'], 'starsolo', '{sample}', '{sample}_starsolo_sce.rds'),
        #        sample = get_sample_names())

rule star_index:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        gtf = config['gtf'],
        fa = config['genome']
    output:
        index_path =  op.join(config['working_dir'] , 'data', 'index', 'star', 'SAindex')
    threads:
        config['nthreads']
    params:
        processing_path = op.join(config['working_dir'], 'data', 'index', 'star/'),
        nthreads = config['nthreads'],
        sjdbOverhang = config['sjdbOverhang'],
        indexNbases = config.get('genomeSAindexNbases', 14)
    log:
        op.join(config['working_dir'], 'logs', 'star_indexing.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'star_indexing.txt')
    shell:
      """
    mkdir -p {params.processing_path}
    # cd {params.processing_path}

    (STAR --runThreadN {params.nthreads} \
     --runMode genomeGenerate \
     --sjdbGTFfile {input.gtf} \
     --genomeDir {params.processing_path} \
     --genomeSAindexNbases {params.indexNbases} \
     --sjdbOverhang {params.sjdbOverhang} \
     --genomeFastaFiles {input.fa} ) 2> {log}
        """


# rule prepare_whitelists:
#     # conda:
#     #     op.join('envs', 'all_in_one.yaml')
#     input:
#         cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
#         cdna = lambda wildcards: get_cdna_by_name(wildcards.sample)
#     output:
#         cb1 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS1.txt'),
#         cb2 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS2.txt'),
#         cb3 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS3.txt')
#     run:
#         sample = wildcards.sample
#         symlink_whitelist(sample)
       
rule starsolo:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        # cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        # cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        standardized_cdna = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cdna.fq.gz"),
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz"),
        index_flag = op.join(config['working_dir'] , 'data', 'index', 'star', 'SAindex'),
        gtf = config['gtf'],
        # cb1 = op.join(config['working_dir'], 'starsolo', "{sample}",  'whitelists', 'BD_CLS1.txt'),
        # cb2 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS2.txt'),
        # cb3 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS3.txt')
    output:
        bam = op.join(config['working_dir'], 'starsolo', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        raw_count_table = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out', 'Gene',
                                  'filtered', 'matrix.mtx'),
        raw_matrix = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out', 'Gene',
                             'raw', 'matrix.mtx')
    threads:
        min(10, config['nthreads'])
    log:
        op.join(config['working_dir'], 'logs', 'starsolo_{sample}.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'star_{sample}.txt')
    params:
        threads = min(10, workflow.cores),
        path = op.join(config['working_dir'], 'starsolo', "{sample}/"),
        index_path = op.join(config['working_dir'] , 'data', 'index', 'star'),
        # num_cells = get_expected_cells_by_name("{sample}"),
        tmp = op.join(config['working_dir'], 'tmp_starsolo_{sample}'),
        maxmem = config['max_mem_mb'] * 1024 * 1024,
        sjdbOverhang = config['sjdbOverhang'],
        soloCellFilter = config['soloCellFilter'],
        soloMultiMappers = config['soloMultiMappers'],
        soloUMIdedup = config.get('soloUMIdedup', '1MM_CR'),
        extraStarSoloArgs = config['extraStarSoloArgs'],
        gene_solo_path = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out',
                            'Gene')
    shell:
        """
   rm -rf {params.tmp} {params.path}/Solo.out
   mkdir -p {params.path} 

   STAR --runThreadN {params.threads} \
        --genomeDir {params.index_path} \
        --readFilesCommand zcat \
        --outFileNamePrefix {params.path} \
        --readFilesIn  {input.standardized_cdna} {input.standardized_cb_umi}  \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 --soloCBlen 27  \
        --soloUMIstart 28 --soloUMIlen 8 \
        --soloUMIdedup {params.soloUMIdedup} \
        --soloBarcodeReadLength 1 \
        --soloCellReadStats Standard \
        --soloCBwhitelist None \
        --soloCellFilter {params.soloCellFilter} \
        --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --sjdbGTFfile {input.gtf} \
        --outTmpDir {params.tmp} \
        --sjdbOverhang {params.sjdbOverhang} \
        --limitBAMsortRAM {params.maxmem} \
        --soloMultiMappers {params.soloMultiMappers} {params.extraStarSoloArgs} 2> {log}

        rm -rf {params.tmp}
        """

# ruleorder: starsolo > symlink_filtered

# rule symlink_filtered:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         bam = op.join(config['working_dir'], 'starsolo', '{sample}', 'Aligned.sortedByCoord.out.bam'),
#         raw_count_table = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out', 'Gene',
#                                   'raw', 'matrix.mtx')
#     output:
#         filtered  = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out',
#                             'Gene', 'filtered', 'matrix.mtx')
#     params:
#         gene_solo_path = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out',
#                             'Gene')
#     threads:
#         1
#     shell:
#         """
#         ## if no cell filtering occurs - happens when heavily downsampling, or when soloCellFilter equals None
#      if [ ! -e {output.filtered} ]; then
#         echo "Caution no cell filtering - symlinking instead"
#         cd {params.gene_solo_path}
#         ln -s {params.gene_solo_path}/raw/*  -t {params.gene_solo_path}/filtered
#      fi
#         """
        
# checkpoint retrieve_genome_sizes:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         fa = config['genome']
#     params:
#         faSize = config['faSize']
#     output:
#         op.join(config['working_dir'], 'data', 'chrom.sizes')
#     shell:
#         """
#         {params.faSize} -detailed -tab {input.fa} > {output}
#         """

# ## TODO keep only if generating coverage tracks
# rule index_bam:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         bam = op.join(config['working_dir'], 'starsolo', '{sample}', 'Aligned.sortedByCoord.out.bam')
#     output:
#         bai = op.join(config['working_dir'], 'starsolo', '{sample}',
#                       'Aligned.sortedByCoord.out.bam.bai')
#     threads: workflow.cores
#     shell:
#         """
#         samtools index -@ {threads} {input.bam}     
#         """


rule generate_sce_starsolo:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        matrix = lambda wildcards: op.join(
            config['working_dir'], 'starsolo', wildcards.sample, 'Solo.out', 'Gene',
            'raw' if config.get('cell_filtering', 'native') == 'emptydrops' else 'filtered',
            'matrix.mtx'),
        bam = op.join(config['working_dir'], 'starsolo', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        script = op.join(config['repo_path'], 'src', 'generate_sce_star.R'),
    output:
        sce = op.join(config['working_dir'], 'starsolo', '{sample}', '{sample}_starsolo_sce.rds')
    params:
        working_dir = config['working_dir'],
        cell_filtering = config.get('cell_filtering', 'native'),
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_star.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_star.txt')
    shell:
        """
        chmod -R ug+rwX $(dirname {input.bam})

        R -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --cell_filtering {params.cell_filtering} \
             --output_fn {output.sce} &> {log}
        """

rule generate_sce_kallisto:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        flag = op.join(config['working_dir'], 'bustools', '{sample}', 'output.mtx'),
        script = op.join(config['repo_path'], 'src', 'generate_sce_kallisto.R'),
    output:
        sce = op.join(config['working_dir'], 'kallisto', '{sample}', '{sample}_kallisto_sce.rds')
    params:
        working_dir = config['working_dir'],
        cell_filtering = config.get('cell_filtering', 'native'),
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_kallisto.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_kallisto.txt')
    shell:
        """
        R -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --cell_filtering {params.cell_filtering} \
             --output_fn {output.sce} &> {log}
        """

rule generate_sce_alevin:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        flag = op.join(config['working_dir'], 'alevin', '{sample}', 'alevin', 'quants_mat.gz'),
        script = op.join(config['repo_path'], 'src', 'generate_sce_alevin.R'),
    output:
        sce = op.join(config['working_dir'], 'alevin', '{sample}', '{sample}_alevin_sce.rds')
    params:
        working_dir = config['working_dir'],
        cell_filtering = config.get('cell_filtering', 'native'),
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_alevin.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_alevin.txt')
    shell:
        """
        R -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --cell_filtering {params.cell_filtering} \
             --output_fn {output.sce} &> {log}
        """

# rule generate_sce_tasseq:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         wta_filtered = op.join(config['working_dir'], 'tasseq', '{sample}', 'Solo.out',
#                                'Gene', 'filtered', 'matrix.mtx'),
#         # gtf = config['gtf'],
#         script = op.join(config['repo_path'], 'src', 'generate_sce.R'),
#         installs = op.join(config['working_dir'], 'logs', 'installs.log')
#     output:
#         sce = op.join(config['working_dir'], 'starsolo', '{sample}', '{sample}_tasseq_sce.rds')
#     params:
#         align_path = op.join(config['working_dir'], 'tasseq'),
#         working_dir = config['working_dir'],
#         sample = lambda wildcards: wildcards.sample,
# #     shell:
#         """
#         R -q --no-save --no-restore --slave \
#              -f {input.script} --args \
#              --sample {wildcards.sample} \
#              --working_dir {params.working_dir} \
#              --output_fn {output.sce}
#         """


rule render_descriptive_report:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        script = op.join(config['repo_path'], 'src', 'generate_descriptive_singlecell_report.Rmd'),
        sces = expand(op.join(config['working_dir'], '{{aligner}}', '{sample}', '{sample}_{{aligner}}_sce.rds'),
                      sample = get_sample_names()),
        counts = ([] if config.get('skip_sampletags', False) else
                  expand(op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_counts.tsv.gz'),
                         sample = get_sample_names()))
    output:
        html = op.join(config['working_dir'], '{aligner}', '{sample}', 'descriptive_report.html')
        # cache = temp(op.join(config['repo_path'], 'process_sce_objects_cache')),
        # cached_files = temp(op.join(config['repo_path'], 'process_sce_objects_files'))
    log:
        op.join(config['working_dir'], 'logs', '{aligner}_{sample}_descriptive_report.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{aligner}_{sample}_descriptive_report.txt')
    params:
        path = op.join(config['working_dir'], '{aligner}', '{sample}'),
        working_dir = op.join(config['working_dir'], '{aligner}'),
        sample = lambda wildcards: wildcards.sample,
    shell:
        """
        cd {params.working_dir}
        mkdir -p {params.path}
        R --vanilla -e 'rmarkdown::render(\"{input.script}\", 
          output_file = \"{output.html}\", 
          params = list(path = \"{params.path}\"))' &> {log}
        """

rule install_rustody:
    conda:
        op.join('envs', 'all_in_one.yaml')
    output:
        op.join('soft', 'Rustody', 'target', '.rustc_info.json')
    log:
        op.join(config['working_dir'], 'logs', 'rustody_install.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'rustody_install.txt')
    shell:
        """
        mkdir -p soft
        cd soft
        curl https://sh.rustup.rs -sSf | sh
        source "$HOME/.cargo/env"
        git clone https://github.com/stela2502/Rustody --depth 1
        git log | head
        cd Rustody
        cargo build --release 2> {log}
        """
        
rule rustody_run:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        rustody = op.join('soft', 'Rustody', 'target', '.rustc_info.json'),
        transcriptome = config['transcriptome'],
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        cb_umi = lambda wildcards: get_cbumi_by_name(wildcards.sample)        
    output:
        flag = op.join(config['working_dir'], 'rustody', '{sample}', 'flag')
    threads:
        workflow.cores
    params:
        whitelist = lambda wildcards: get_barcode_whitelist_by_name(wildcards.sample),
        species = lambda wildcards: get_species_by_name(wildcards.sample),
        rustody_path = op.join('soft', 'Rustody', 'target', 'release')
    log:
        op.join(config['working_dir'], 'logs', 'rustody_{sample}.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'rustody_{sample}.txt')
    shell:
        """
        source "$HOME/.cargo/env"
        if [ {params.whitelist} == '384x3' ]; then
           wl='v2.384'
        elif [ {params.whitelist} == '96x3' ]; then
           wl='v2.96'
        else
           wl='error_unknown_whitelist_spec'
        fi

        export PATH="{params.rustody_path}:"$PATH
        quantify_rhapsody_multi \
           --version "$wl" \
           --specie {params.species} \
           --reads {input.cb_umi} \
           --outpath ~/Rustody_pdgfra \
           --num-threads {threads} \
           --file {input.cdna} \
           --expression {input.transcriptome} \
           --min-umi 1 &> {log}

        touch {output.flag}
        """

rule standardize_cb_umis_cutadapt:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        cb_umi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
        cdna = lambda wildcards: get_cdna_by_name(wildcards.sample)
    output:
        cdna = op.join(config['working_dir'], 'data', 'fastq',
                       "{sample}_standardized_cdna.fq.gz"),
        temp_cb_umi = temp(op.join(config['working_dir'], 'data', 'fastq',
                                      "{sample}_temp_cb_umi.fq.gz")),        
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq',
                                      "{sample}_standardized_cb_umi.fq.gz")
    params:
        path = op.join(config['working_dir'], 'data', 'fastq')
    threads:
        workflow.cores
    log:
        op.join(config['working_dir'], 'logs', 'standardize_cb_umis_{sample}.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'standardize_cb_umis_{sample}.txt')
    shell:
        """
        mkdir -p {params.path}
        cutadapt -g "NNNNNNNNNGTGANNNNNNNNNGACANNNNNNNNNNNNNNNNN;min_overlap=43;noindels" \
           --action=crop \
           --discard-untrimmed \
           --cores {threads} \
           -e 0 \
           --pair-filter=both \
           -o {output.temp_cb_umi} \
           -p {output.cdna} \
           {input.cb_umi} {input.cdna} &> {log}

        pigz --decompress {output.temp_cb_umi} -p {threads} --stdout | \
            cut -c1-9,14-22,27- | pigz -p {threads} > {output.standardized_cb_umi}
        """

rule kallisto_index:
    conda:
        op.join('envs', 'kallisto.yaml')
    input:
        transcriptome = config['transcriptome'],
    params:
        index_name = 'kallisto.index',
        output_dir= op.join(config['working_dir'], 'data', 'index', 'kallisto')        
    output:
        index_name = op.join(config['working_dir'], 'data', 'index', 'kallisto', 'kallisto.index')
    threads:
        workflow.cores
    log:
        op.join(config['working_dir'], 'logs', 'kallisto_index.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'kallisto_index.txt')
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        
        kallisto index --threads {threads} -i {params.index_name} {input.transcriptome} &> {log}
        """

rule kallisto_bus:
    conda:
        op.join('envs', 'kallisto.yaml')
    input:
        transcriptome = config['transcriptome'],
        # transcriptome = op.join(config['working_dir'], 'data', 'index', 'salmon', 'transcriptome.fa'),        
        # cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        standardized_cdna = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cdna.fq.gz"),
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz"),
        kallisto_index = op.join(config['working_dir'], 'data', 'index', 'kallisto', 'kallisto.index'),
    output:
        matrix_ec = op.join(config['working_dir'], 'kallisto', '{sample}', 'matrix.ec'),
        transcripts = op.join(config['working_dir'], 'kallisto', '{sample}', 'transcripts.txt'),
        bus = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.bus'),
        tmp = temp(op.join(config['working_dir'], 'kallisto', '{sample}', 'transcripts.txt.wrong'))
    params:
        output_dir = op.join(config['working_dir'], 'kallisto', '{sample}'),
        gtf_style = config['gtf_origin']
    threads: workflow.cores
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_kallisto_bus.txt')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_kallisto_bus.log')
    shell:
        """
        kallisto bus --index {input.kallisto_index} \
            --output-dir {params.output_dir} \
            -x '0,0,27:0,27,35:1,0,0' \
            -t {threads} \
            {input.standardized_cb_umi} {input.standardized_cdna} &> {log}

         if [[ {params.gtf_style} == 'ensembl' ]]
         then
            echo "Ensembl GTF, nothing to do"
            touch {output.tmp}   
         elif [[ {params.gtf_style} == 'gencode' ]]
         then
            echo "Gencode GTF, standardizing transcripts"
            sed  's/|/ /g' {output.transcripts} | cut -f1 -d" " > {output.tmp}
            cp {output.tmp} {output.transcripts}
         
         else
           echo "gtf_origin is misspecified within the config file"
         fi  
        """

        
rule bustools_sort:
    conda:
        op.join('envs', 'kallisto.yaml')
    input:
        bus    = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.bus')
    output:
        sorted_bus = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.sorted.bus')
    threads: workflow.cores
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_bustools_sort.txt')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_bustools_sort.log')
    shell:
        """
        bustools sort -t {threads} -o {output.sorted_bus} {input.bus} &> {log}
        """

rule derive_kallisto_observed_whitelist:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        sorted_bus = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.sorted.bus'),
        cb_umi_fq = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz"),
        cb1 = lambda wildcards: op.join(config['repo_path'], 'data', 'whitelist_' + get_barcode_whitelist_by_name(wildcards.sample), 'BD_CLS1.txt'),
        cb2 = lambda wildcards: op.join(config['repo_path'], 'data', 'whitelist_' + get_barcode_whitelist_by_name(wildcards.sample), 'BD_CLS2.txt'),
        cb3 = lambda wildcards: op.join(config['repo_path'], 'data', 'whitelist_' + get_barcode_whitelist_by_name(wildcards.sample), 'BD_CLS3.txt')
    output:
        observed_whitelist = op.join(config['working_dir'], 'kallisto', '{sample}', 'observed_whitelist.txt')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_derive_observed_whitelist.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_derive_observed_whitelist.txt')
    shell:
        """
        set -euo pipefail
        python - <<'PY' > {output.observed_whitelist}
import gzip

with open("{input.cb1}") as fh:
    cb1 = {line.strip() for line in fh if line.strip()}
with open("{input.cb2}") as fh:
    cb2 = {line.strip() for line in fh if line.strip()}
with open("{input.cb3}") as fh:
    cb3 = {line.strip() for line in fh if line.strip()}

observed_valid_barcodes = set()
with gzip.open("{input.cb_umi_fq}", "rt") as fh:
    for line_index, line in enumerate(fh):
        if line_index % 4 != 1:
            continue
        sequence = line.strip()
        if len(sequence) < 27:
            continue
        barcode = sequence[:27]
        barcode_1 = barcode[:9]
        barcode_2 = barcode[9:18]
        barcode_3 = barcode[18:27]
        if barcode_1 in cb1 and barcode_2 in cb2 and barcode_3 in cb3:
            observed_valid_barcodes.add(barcode)

for barcode in sorted(observed_valid_barcodes):
    print(barcode)
PY
        wc -l {output.observed_whitelist} > {log}
        """

rule bustools_correct:
    conda:
        op.join('envs', 'kallisto.yaml')
    input:
        bus = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.sorted.bus'),
        whitelist = (
            op.join(config['working_dir'], 'simulate', 'cell_barcodes.txt')
            if config.get('use_simulated', False)
            else op.join(config['working_dir'], 'kallisto', '{sample}', 'observed_whitelist.txt')
        )
    output:
        corrected_bus = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.corrected.bus')
    threads: 1
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_bustools_correct.txt')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_bustools_correct.log')
    shell:
        """
        bustools correct -w {input.whitelist} -o {output.corrected_bus} {input.bus} &> {log}
        """

rule bustools_count:
    conda:
        op.join('envs', 'kallisto.yaml')
    input:
        txp2gene   = op.join(config['working_dir'], 'data', 'index', 'salmon', 'txp2gene'),
        matrix_ec  = op.join(config['working_dir'], 'kallisto', '{sample}', 'matrix.ec'),
        transcripts = op.join(config['working_dir'], 'kallisto', '{sample}', 'transcripts.txt'),
        bus        = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.corrected.bus')
    output:
        op.join(config['working_dir'], 'bustools', '{sample}', 'output.mtx')
    params:
        output_dir = op.join(config['working_dir'], 'bustools', '{sample}/')
    threads: 1
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_bustools_count.txt')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_bustools_count.log')
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        bustools count \
           -o {params.output_dir} \
           -g {input.txp2gene} \
           -e {input.matrix_ec} \
           -t {input.transcripts} \
           --genecounts \
           {input.bus} &> {log}
        """
        

        
## from https://github.com/imallona/rock_roi_paper/blob/imallona/03_leukemia/02_sampletags_again.sh

for sample in get_sample_names():
    species = get_species_by_name(name = sample)
    rule:
        name:
            f"{species}_star_index_sampletags"
        conda:
            op.join('envs', 'all_in_one.yaml')
        input:
            fa = op.join('data', 'sampletags', species + '_sampletags.fa')
        output:
            op.join(config['working_dir'] , 'data', species + '_index', 'sampletags', 'SAindex')
        threads:
            workflow.cores
        params:
            output_dir = op.join(config['working_dir'], 'data', species + '_index', 'sampletags')
        log:
            op.join(config['working_dir'], 'logs', species + '_sampletags_index.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', species + '_sampletags__index.txt')
        shell:
            """
            STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeSAindexNbases 2 \
            --genomeDir {params.output_dir} \
            --genomeFastaFiles {input.fa} &> {log}
            """
        
rule extract_unmapped_startsolo_wta_tagged_fastqs:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        bam = op.join(config['working_dir'], 'starsolo', '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        temp(op.join(config['working_dir'], 'sampletags', '{sample}_unmapped_tagged.fq.gz'))
    threads:
        min(10, workflow.cores)
    log:
        op.join(config['working_dir'], 'logs', '{sample}_unmapped_tagged.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_unaligned_tagged.txt')
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view -@ {threads} -f 4 -d CB {input.bam} | \
            awk -F '\\t' '{{
                cb=""; ub="";
                for (i=12; i<=NF; i++) {{
                    if (substr($i,1,5)=="CB:Z:") cb=substr($i,6);
                    else if (substr($i,1,5)=="UB:Z:") ub=substr($i,6);
                }}
                if (cb!="" && ub!="") printf "@%s__%s\\n%s\\n+\\n%s\\n", cb, ub, $10, $11
            }}' | pigz -p {threads} -c > {output} 2> {log}
        """
    

rule extract_sampletagslooking_fastqs:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        fq = op.join(config['working_dir'], 'sampletags', '{sample}_unmapped_tagged.fq.gz')
    output:
        fq = temp(op.join(config['working_dir'], 'sampletags', '{sample}_sampletag_tagged.fq.gz'))
    threads:
        workflow.cores
    log:
        op.join(config['working_dir'], 'logs', '{sample}_sampletags_cutadapt.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_sampletags_cutadapt.txt')
    shell:
        """
        cutadapt -g ^GTTGTCAAGATGCTACCGTTCAGAG {input.fq} \
          -j {threads} --action=retain --discard-untrimmed \
          -o {output.fq} &> {log}
        """
        
rule align_star_sampletags:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        fastq = op.join(config['working_dir'], 'sampletags', '{sample}_sampletag_tagged.fq.gz'),
        idx_flag = lambda wildcards: op.join(config['working_dir'], 'data',
                       get_species_by_name(wildcards.sample) + '_index', 'sampletags', 'SAindex')
    output:
        bam = op.join(config['working_dir'], 'sampletags', '{sample}', 'Aligned.out.bam')
    params:
        output_dir = op.join(config['working_dir'], 'sampletags', '{sample}/'),
        sampletags_genome_dir = lambda wildcards: op.join(config['working_dir'], 'data',
                                    get_species_by_name(wildcards.sample) + '_index', 'sampletags'),
        tmp = op.join(config['working_dir'], 'sampletags', 'tmp_starsolo_{sample}'),
    log:
        op.join(config['working_dir'], 'logs', '{sample}_align_sampletags_star.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_align_sampletags_star.txt')
    threads:
        min(10, workflow.cores)
    shell:
        """
        rm -rf {params.tmp} {params.output_dir}
        mkdir -p {params.output_dir} 
        
        STAR --runThreadN {threads} \
          --genomeDir {params.sampletags_genome_dir} \
          --outTmpDir {params.tmp} \
          --readFilesCommand zcat \
          --outFileNamePrefix {params.output_dir} \
          --readFilesIn {input.fastq}  \
          --outSAMtype BAM Unsorted \
          --outFilterMismatchNmax 5 \
          --scoreInsOpen  -8 \
          --scoreDelBase -8 \
          --alignIntronMax 1 &> {log}

        rm -rf {params.tmp}
        """

rule count_sampletags:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        bam = op.join(config['working_dir'], 'sampletags', '{sample}', 'Aligned.out.bam')
    output:
        counts = op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_counts.tsv.gz')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_sampletag_counting.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_sampletag_counting.txt')
    threads:
        min(10, workflow.cores)
    shell:
        """
        ## this only reports a table with as many rows as `cb,umi,sampletag,cigar` alignments. Not summarized
        ##  at all 
        samtools view -@ {threads} {input.bam} | \
          cut -f1,3,6 | sed 's/__/\t/g' | pigz -p {threads} -c > {output.counts}
        """

rule render_sampletag_report:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        counts   = op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_counts.tsv.gz'),
        script   = op.join(config['repo_path'], 'src', 'generate_sampletag_report.Rmd'),
    output:
        html = op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_report.html')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_sampletag_report.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_sampletag_report.txt')
    params:
        assignments = lambda wildcards: (
            op.join(config['working_dir'], 'simulate', 'sampletag_assignments.txt')
            if config.get('use_simulated', False) else ''
        )
    shell:
        """
        R --vanilla -e \
          'rmarkdown::render("{input.script}",
            output_file = "{output.html}",
            params = list(
              counts_path      = "{input.counts}",
              assignments_path = "{params.assignments}"))' &> {log}
        """

rule deversion_transcriptome:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        transcriptome = config['transcriptome']
    output:
        deversioned_fasta = op.join(config['working_dir'], 'data', 'index', 'salmon', 'transcriptome.fa')
    params:
        index_path = op.join(config['working_dir'], 'data', 'index', 'salmon'),
        gtf_style = config['gtf_origin']
    threads: workflow.cores    
    log:
        op.join(config['working_dir'], 'logs', 'deversion_transcriptome.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'deversion_transcriptome.txt')
    shell:
        """
        mkdir -p {params.index_path}

         if [[ {params.gtf_style} == 'ensembl' ]]
         then
            echo "Ensembl GTF"

            # no transcript versions
            #  e.g. ENST4654.1, the .1 needs to go because the matching GTF doesn't have it
            zcat {input.transcriptome} | awk '/^>/ {{split($1, a, "."); $1=a[1]}} {{print}}' > {output.deversioned_fasta}

         elif [[ {params.gtf_style} == 'gencode' ]]
         then
            echo "Gencode GTF"
            ## yes versions, nbo pipes
            zcat {input.transcriptome} | sed 's/|/ /g' > {output.deversioned_fasta}

         else
           echo "gtf_origin is misspecified within the config file"
         fi      
        """
    
rule salmon_index:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        transcriptome = op.join(config['working_dir'], 'data', 'index', 'salmon', 'transcriptome.fa')
    output:
        index_flag = op.join(config['working_dir'], 'data', 'index', 'salmon', 'seq.bin')
    params:
        index_path = op.join(config['working_dir'], 'data', 'index', 'salmon'),
        gtf_style = config['gtf_origin']
    threads: workflow.cores    
    log:
        op.join(config['working_dir'], 'logs', 'salmon_index.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'salmon_index.txt')
    shell:
        """
        mkdir -p {params.index_path}

        salmon index -t {input.transcriptome} -i {params.index_path} -p {threads} &> {log}        
        """
         
rule get_txp2gene:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        genes_gtf = config['gtf']
    output:
        op.join(config['working_dir'], 'data', 'index', 'salmon', 'txp2gene')
    params:
        gtf_style = config['gtf_origin']
    log:
        op.join(config['working_dir'], 'logs', 'get_txp2gene.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'get_txp2gene.txt')
    shell:
         """
         if [[ {params.gtf_style} == 'ensembl' ]]
         then
            echo "Ensembl GTF"
            grep transcript {input.genes_gtf} | \
                 awk '{{print $14,$10}}' | sed -e 's|"||g' -e 's|;||g' | uniq > {output}
         elif [[ {params.gtf_style} == 'gencode' ]]
         then
            echo "Gencode GTF"
            grep transcript {input.genes_gtf} | \
                 awk '{{print $12,$10}}' | sed -e 's|"||g' -e 's|;||g' | uniq > {output}
         else
           echo "gtf_origin is misspecified within the config file"
         fi         
         """
         
rule alevin_align:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        index_flag = op.join(config['working_dir'], 'data', 'index', 'salmon', 'seq.bin'),
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq',
                                      "{sample}_standardized_cb_umi.fq.gz"),
        # cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        standardized_cdna = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cdna.fq.gz"),
        t2g = op.join(config['working_dir'], 'data', 'index', 'salmon', 'txp2gene')
    params:
        index_path =  op.join(config['working_dir'], 'data', 'index', 'salmon'),
        output_dir = op.join(config['working_dir'], 'alevin', '{sample}')
    output:
        flag = op.join(config['working_dir'], 'alevin', '{sample}', 'alevin', 'quants_mat.gz')
    threads:
        workflow.cores      
    log:
        op.join(config['working_dir'], 'logs', 'alevin_{sample}_align.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'alevin_{sample}_align.txt')
    shell:
        """
        salmon alevin -i {params.index_path} \
           -l ISR \
           -1 {input.standardized_cb_umi} \
           -2 {input.standardized_cdna} \
           --bc-geometry '1[1-27]' \
           --read-geometry '2[1-end]' \
           --umi-geometry '1[28-35]' \
           -o {params.output_dir} \
           -p {threads} \
           --tgMap {input.t2g} &> {log}
        """


## SevenBridges / BD Rhapsody official pipeline (optional)
##
## Configure in the config yaml, either globally at the top level or per-sample
## inside uses:. Per-sample values override the global ones.
##
## Run mode (sbg_cwl is set):
##   Executes the BD Rhapsody CWL workflow via cwl-runner --singularity.
##   cwl-runner pulls the BD Docker image listed inside the CWL file from
##   DockerHub and converts it to a Singularity image automatically; no
##   manual docker/singularity pull is needed.
##   The CWL input.yml is generated entirely from config values; users never
##   edit a separate yml file.
##
##   Reference archive handling:
##   - Simulated data (use_simulated: true): assembled from the simulated STAR
##     index and GTF following the BD Rhapsody format:
##       BD_Rhapsody_Reference_Files/
##         star_index/   [STAR genomeGenerate output]
##         *.gtf
##   - Real data, sbg_reference_url set: downloaded from the BD public S3 bucket
##       http://bd-rhapsody-public.s3-website-us-east-1.amazonaws.com/Rhapsody-WTA/
##   - Real data, sbg_reference_archive set: path to a pre-built local archive
##   - Real data, neither set: built automatically from the STAR index and GTF
##
## generate_sce_sbg decodes numeric BD barcode indices to 27-bp sequences
## via scripts/index2barcode.R.

## Simulated reference: built from the simulated STAR index + GTF.
## Required format documented at:
## https://bd-rhapsody-bioinfo-docs.genomics.bd.com/setup/input/reference_files.html
if config.get('use_simulated') and any(get_sbg_cwl_by_name(s) for s in get_sample_names()):
    rule build_sbg_reference_simulated:
        conda:
            op.join('envs', 'all_in_one.yaml')
        input:
            star_flag = op.join(config['working_dir'], 'data', 'index', 'star', 'SAindex'),
            gtf = config['gtf']
        output:
            archive = op.join(config['working_dir'], 'sbg_reference',
                              'rhapsody_reference_simulated.tar.gz')
        params:
            star_index_dir = op.join(config['working_dir'], 'data', 'index', 'star'),
            staging = op.join(config['working_dir'], 'sbg_reference', 'staging')
        log:
            op.join(config['working_dir'], 'logs', 'build_sbg_reference_simulated.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', 'build_sbg_reference_simulated.txt')
        shell:
            """
            rm -rf {params.staging}
            mkdir -p {params.staging}/BD_Rhapsody_Reference_Files/star_index
            cp {params.star_index_dir}/* {params.staging}/BD_Rhapsody_Reference_Files/star_index/
            cp {input.gtf} {params.staging}/BD_Rhapsody_Reference_Files/
            tar -czf {output.archive} -C {params.staging} BD_Rhapsody_Reference_Files 2> {log}
            rm -rf {params.staging}
            """


## Real-data reference: three modes (mutually exclusive, evaluated at parse time):
##   1. sbg_reference_url set   → download_sbg_reference (curl from URL)
##   2. sbg_reference_archive set → use pre-built archive, no Snakemake tracking
##   3. neither set             → build_sbg_reference_real (package STAR index + GTF)
_any_sbg_ref_url = (
    config.get('sbg_reference_url') or
    any(_sbg_uses(s, 'sbg_reference_url') for s in get_sample_names())
)
_global_sbg_ref_url = config.get('sbg_reference_url') or next(
    (_sbg_uses(s, 'sbg_reference_url') for s in get_sample_names()
     if _sbg_uses(s, 'sbg_reference_url')),
    ''
)
_any_sbg_ref_archive = any(
    get_sbg_reference_by_name(s) for s in get_sample_names()
)

if _any_sbg_ref_url and not config.get('use_simulated'):
    rule download_sbg_reference:
        output:
            archive = op.join(config['working_dir'], 'sbg_reference',
                              'rhapsody_reference.tar.gz')
        params:
            url = _global_sbg_ref_url
        log:
            op.join(config['working_dir'], 'logs', 'download_sbg_reference.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', 'download_sbg_reference.txt')
        shell:
            """
            mkdir -p $(dirname {output.archive})
            curl -fsSL -o {output.archive} "{params.url}" &> {log}
            """

## Build the SBG reference locally from the existing STAR index and GTF.
## Same layout as build_sbg_reference_simulated; triggered when sbg_cwl is
## configured but neither sbg_reference_url nor sbg_reference_archive is set.
if (not config.get('use_simulated') and _has_sbg
        and not _any_sbg_ref_url and not _any_sbg_ref_archive):
    rule build_sbg_reference_real:
        conda:
            op.join('envs', 'all_in_one.yaml')
        input:
            star_flag = op.join(config['working_dir'], 'data', 'index', 'star', 'SAindex'),
            gtf = config['gtf']
        output:
            archive = op.join(config['working_dir'], 'sbg_reference',
                              'rhapsody_reference.tar.gz')
        params:
            star_index_dir = op.join(config['working_dir'], 'data', 'index', 'star'),
            staging = op.join(config['working_dir'], 'sbg_reference', 'staging')
        log:
            op.join(config['working_dir'], 'logs', 'build_sbg_reference_real.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', 'build_sbg_reference_real.txt')
        shell:
            """
            rm -rf {params.staging}
            mkdir -p {params.staging}/BD_Rhapsody_Reference_Files/star_index
            cp {params.star_index_dir}/* {params.staging}/BD_Rhapsody_Reference_Files/star_index/
            cp {input.gtf} {params.staging}/BD_Rhapsody_Reference_Files/
            tar -czf {output.archive} -C {params.staging} BD_Rhapsody_Reference_Files 2> {log}
            rm -rf {params.staging}
            """


## run_sbg_cwl fires for every sample when 'sbg' is in the aligner list.
## At parse time we decide which reference file to track based on the mode.
## caution it provides the unfiltered counts as filtered ones if the filtering fails
for _sbg_sample in _sbg_samples:

    _sbg_ref_url = (
        _sbg_uses(_sbg_sample, 'sbg_reference_url') or config.get('sbg_reference_url')
    )

    if config.get('use_simulated'):
        _ref_path = op.join(config['working_dir'], 'sbg_reference',
                            'rhapsody_reference_simulated.tar.gz')
        _ref_input = [_ref_path]
    elif _sbg_ref_url:
        _ref_path = op.join(config['working_dir'], 'sbg_reference',
                            'rhapsody_reference.tar.gz')
        _ref_input = [_ref_path]
    elif get_sbg_reference_by_name(_sbg_sample):
        _ref_path = get_sbg_reference_by_name(_sbg_sample)
        _ref_input = []  # pre-existing archive; not tracked by Snakemake
    else:
        ## auto-build from local STAR index + GTF (build_sbg_reference_real)
        _ref_path = op.join(config['working_dir'], 'sbg_reference',
                            'rhapsody_reference.tar.gz')
        _ref_input = [_ref_path]

    rule:
        name: f"run_sbg_cwl_{_sbg_sample}"
        conda:
            op.join('envs', 'sbg.yaml')
        input:
            r1 = get_cbumi_by_name(_sbg_sample),
            r2 = get_cdna_by_name(_sbg_sample),
            ref = _ref_input
        output:
            matrix_unfiltered = op.join(config['working_dir'], 'sbg', _sbg_sample,
                                        'unfiltered_MEX_output', 'matrix.mtx.gz'),
            matrix_filtered = op.join(config['working_dir'], 'sbg', _sbg_sample,
                                      'filtered_MEX_output', 'matrix.mtx.gz')
        params:
            cwl = get_sbg_cwl_by_name(_sbg_sample),
            outdir = op.join(config['working_dir'], 'sbg', _sbg_sample),
            sample_tags_version = get_sbg_sample_tags_version_by_name(_sbg_sample),
            ref_path = _ref_path
        log:
            op.join(config['working_dir'], 'logs', f'sbg_cwl_{_sbg_sample}.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', f'sbg_cwl_{_sbg_sample}.txt')
        shell:
            """
            mkdir -p {params.outdir}

            ## Generate the CWL input.yml entirely from config values.
            INPUT_YML={params.outdir}/input.yml
            cat > "$INPUT_YML" << 'ENDOFYML'
Reads:
  - class: File
    location: "PLACEHOLDER_R1"
  - class: File
    location: "PLACEHOLDER_R2"

Reference_Archive:
    class: File
    location: "PLACEHOLDER_REF"

Sample_Tags_Version: PLACEHOLDER_STV

## Force short-read (STAR) mode. Without this, v3.0 auto-detects and may
## switch to bwa-mem2 for long reads, which fails when the reference
## archive has no bwa-mem2 index (only a STAR index).
Long_Reads: false
ENDOFYML

            sed -i "s|PLACEHOLDER_R1|$(realpath {input.r1})|" "$INPUT_YML"
            sed -i "s|PLACEHOLDER_R2|$(realpath {input.r2})|" "$INPUT_YML"
            sed -i "s|PLACEHOLDER_REF|$(realpath {params.ref_path})|" "$INPUT_YML"
            sed -i "s|PLACEHOLDER_STV|{params.sample_tags_version}|" "$INPUT_YML"

            cwl-runner --singularity \
                --outdir {params.outdir} \
                {params.cwl} "$INPUT_YML" &> {log}

            ## cwl-runner places MEX outputs as zips in --outdir.
            ## Unzip both filtered and unfiltered into their respective subdirectories.
            UNFILTERED_ZIP=$(ls {params.outdir}/*_RSEC_MolsPerCell_Unfiltered_MEX.zip 2>/dev/null | head -1)
            if [ -z "$UNFILTERED_ZIP" ] || [ ! -f "$UNFILTERED_ZIP" ]; then
                echo "ERROR: expected unfiltered MEX zip not found in {params.outdir}" >> {log}
                ls -lR {params.outdir} >> {log} 2>&1
                exit 1
            fi
            mkdir -p {params.outdir}/unfiltered_MEX_output
            unzip -o "$UNFILTERED_ZIP" -d {params.outdir}/unfiltered_MEX_output >> {log} 2>&1

            FILTERED_ZIP=$(ls {params.outdir}/*_RSEC_MolsPerCell_MEX.zip 2>/dev/null | grep -v Unfiltered | head -1)
            if [ -z "$FILTERED_ZIP" ] || [ ! -f "$FILTERED_ZIP" ]; then
                echo "WARNING: expected filtered MEX zip not found in {params.outdir}, using unfiltered instead" >> {log}
                
                ## fallback to link the unfiltered as filtered
                FILTERED_ZIP="${{UNFILTERED_ZIP/Unfiltered_/Filtered_}}"
                ln -sf "$UNFILTERED_ZIP" "$FILTERED_ZIP"
            
            fi
            mkdir -p {params.outdir}/filtered_MEX_output
            unzip -o "$FILTERED_ZIP" -d {params.outdir}/filtered_MEX_output >> {log} 2>&1
            """


if _has_sbg:
    _sbg_mex_subdir = ('unfiltered_MEX_output'
                       if config.get('cell_filtering', 'native') == 'emptydrops'
                       else 'filtered_MEX_output')

    rule generate_sce_sbg:
        conda:
            op.join('envs', 'all_in_one.yaml')
        input:
            mex_flag = lambda wildcards: op.join(
                config['working_dir'], 'sbg', wildcards.sample,
                _sbg_mex_subdir, 'matrix.mtx.gz'),
            script = op.join(config['repo_path'], 'src', 'generate_sce_sbg.R'),
            index2barcode = op.join(config['repo_path'], 'src', 'index2barcode.R')
        output:
            sce = op.join(config['working_dir'], 'sbg', '{sample}', '{sample}_sbg_sce.rds')
        params:
            mex_dir = lambda wildcards: op.join(
                config['working_dir'], 'sbg', wildcards.sample, _sbg_mex_subdir),
            cell_filtering = config.get('cell_filtering', 'native'),
            bead_version = lambda wildcards: get_sbg_bead_version_by_name(wildcards.sample),
            whitelist_dir = lambda wildcards: (
                op.join(config['repo_path'], 'data',
                        'whitelist_' + get_barcode_whitelist_by_name(wildcards.sample))
                if get_sbg_bead_version_by_name(wildcards.sample) == 'EnhV2' else ''
            ),
            features_map = lambda wildcards: op.join(
                config['working_dir'], 'starsolo', wildcards.sample,
                'Solo.out', 'Gene', 'filtered', 'features.tsv'),
        log:
            op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_sbg.log')
        benchmark:
            op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_sbg.txt')
        shell:
            """
            R -q --no-save --no-restore --slave \
                 -f {input.script} --args \
                 --sample {wildcards.sample} \
                 --mex_dir {params.mex_dir} \
                 --output_fn {output.sce} \
                 --bead_version {params.bead_version} \
                 --index2barcode_script {input.index2barcode} \
                 --cell_filtering {params.cell_filtering} \
                 $([ -n "{params.whitelist_dir}" ] && echo "--whitelist_dir {params.whitelist_dir}") \
                 $([ -f "{params.features_map}" ] && echo "--features_map {params.features_map}") \
                 &> {log}
            """


## Paper-figure reports

rule render_comparison_report:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        sces = lambda wildcards: (
            expand(
                op.join(config['working_dir'], '{aligner}', wildcards.sample,
                        wildcards.sample + '_{aligner}_sce.rds'),
                aligner = get_aligners()
            ) + (
                [op.join(config['working_dir'], 'sbg', wildcards.sample,
                         wildcards.sample + '_sbg_sce.rds')]
                if _has_sbg else []
            )
        ),
        barcodes = (
            op.join(config['working_dir'], 'simulate', 'cell_barcodes.txt')
            if config.get('use_simulated', False) else []
        ),
        true_mex = (
            op.join(config['working_dir'], 'simulate', 'true_mex', 'matrix.mtx.gz')
            if config.get('use_simulated', False) else []
        ),
        doc = op.join(config['repo_path'], 'docs', '02_comparison.Rmd'),
    output:
        html = op.join(config['working_dir'], '{sample}_comparison.html')
    params:
        working_dir = config['working_dir'],
        sample = lambda wildcards: wildcards.sample,
        aligners = ','.join(get_aligners()),
        has_sbg = 'true' if _has_sbg else 'false',
        n_expected_cells = (config.get('sim_n_cells', 0)
                            if config.get('use_simulated', False) else 0),
        n_umis_per_cell = (config.get('sim_n_umis', 0)
                           if config.get('use_simulated', False) else 0),
        barcodes_file = (op.join(config['working_dir'], 'simulate', 'cell_barcodes.txt')
                         if config.get('use_simulated', False) else ''),
        true_mex_dir = (op.join(config['working_dir'], 'simulate', 'true_mex')
                        if config.get('use_simulated', False) else ''),
    log:
        op.join(config['working_dir'], 'logs', '{sample}_comparison_report.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_comparison_report.txt')
    shell:
        """
        R --vanilla -e '
          rmarkdown::render(
            "{input.doc}",
            output_file  = "{output.html}",
            params = list(
              working_dir      = "{params.working_dir}",
              sample           = "{params.sample}",
              aligners         = "{params.aligners}",
              has_sbg          = "{params.has_sbg}",
              n_expected_cells = {params.n_expected_cells},
              n_umis_per_cell  = {params.n_umis_per_cell},
              barcodes_file    = "{params.barcodes_file}",
              true_mex_dir     = "{params.true_mex_dir}"))' &> {log}
        """


rule render_benchmarks_report:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        ## depend on descriptive reports so all benchmark files have been written
        descriptive = expand(
            op.join(config['working_dir'], '{aligner}', '{sample}', 'descriptive_report.html'),
            aligner = get_aligners(),
            sample  = get_sample_names()
        ),
        sbg_sce = _sbg_sce_targets,
        doc = op.join(config['repo_path'], 'docs', '03_benchmarks.Rmd'),
    output:
        html = op.join(config['working_dir'], 'benchmarks_report.html')
    params:
        working_dir = config['working_dir'],
        n_expected_cells = (config.get('sim_n_cells', 0)
                            if config.get('use_simulated', False) else 0),
    log:
        op.join(config['working_dir'], 'logs', 'benchmarks_report.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'benchmarks_report.txt')
    shell:
        """
        R --vanilla -e '
          rmarkdown::render(
            "{input.doc}",
            output_file  = "{output.html}",
            params = list(
              working_dir      = "{params.working_dir}",
              n_expected_cells = {params.n_expected_cells}))' &> {log}
        """


# # https://github.com/s-shichino1989/TASSeq_EnhancedBeads/blob/e48fd2c2fd5a23d622f03e206b8fbe87772fd57f/shell_scripts/Rhapsody_STARsolo.sh#L18
# rule starsolo_tasseq_style:
#     conda:
#         op.join('envs', 'all_in_one.yaml')
#     input:
#         cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
#         cbumi = lambda wildcards: get_cbumi_by_name(wildcards.sample),
#         index_flag = op.join(config['working_dir'] , 'data', 'index', 'SAindex'),
#         gtf = config['gtf'],
#         cb1 = op.join(config['working_dir'], 'starsolo', "{sample}",  'whitelists', 'BD_CLS1.txt'),
#         cb2 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS2.txt'),
#         cb3 = op.join(config['working_dir'], 'starsolo', "{sample}", 'whitelists', 'BD_CLS3.txt')
#     output:
#         bam = op.join(config['working_dir'], 'tasseq', '{sample}', 'Aligned.sortedByCoord.out.bam'),
#         filtered_barcodes = op.join(config['working_dir'], 'tasseq', '{sample}', 'Solo.out', 'Gene',
#                                     'filtered', 'barcodes.tsv'),
#         filtered_counts = op.join(config['working_dir'], 'tasseq', '{sample}', 'Solo.out', 'Gene',
#                                   'filtered', 'matrix.mtx')
#     threads: workflow.cores
#     params:
#         threads = min(10, workflow.cores),
#         path = op.join(config['working_dir'], 'tasseq', "{sample}/"),
#         index_path = op.join(config['working_dir'] , 'data', 'index'),
#         # num_cells = get_expected_cells_by_name("{sample}"),
#         tmp = op.join(config['working_dir'], 'tmp_tasseq_{sample}'),
#         maxmem = config['max_mem_mb'] * 1024 * 1024,
#         sjdbOverhang = config['sjdbOverhang'],
#         soloCellFilter = config['soloCellFilter']
#     shell:
#         """
#         rm -rf {params.tmp}
#         mkdir -p {params.path} 

#         STAR --runThreadN {params.threads} \
#           --genomeDir {params.index_path} \
#         --readFilesIn {input.cdna} {input.cbumi} \
#         --outFileNamePrefix {params.path} \
#         --readFilesCommand zcat \
#         --clipAdapterType CellRanger4 \
#         --outSAMtype BAM SortedByCoordinate \
#         --outBAMsortingThreadN {params.threads} \
#         --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
#         --outSAMunmapped Within \
#         --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
#         --outFilterMultimapScoreRange 0 --seedSearchStartLmax 30 \
#         --soloCellFilter {params.soloCellFilter} \
#         --soloUMIdedup Exact \
#         --soloMultiMappers Rescue \
#         --soloFeatures Gene GeneFull \
#         --soloAdapterSequence NNNNNNNNNGTGANNNNNNNNNGACA \
#         --soloCBmatchWLtype EditDist_2 \
#         --soloCBwhitelist {input.cb1} {input.cb2} {input.cb3} \
#         --soloType CB_UMI_Complex \
#         --soloUMIlen 8 \
#         --soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 \
#         --soloUMIposition 3_10_3_17 
#     rm -rf {params.tmp}
#         """
