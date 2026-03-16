#!/usr/bin/env snakemake -s
##
## Snakefile to process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

import os.path as op
import os

configfile: "config.yaml"

## to ease whitelists symlinking — must be absolute before any include uses it
if not op.isabs(config['repo_path']):
    config['repo_path'] = op.join(workflow.basedir, config['repo_path'])

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

## kallisto and bustools from bioconda are not reliable, so we compile them
kallisto = op.join(config['working_dir'], 'software', 'kallisto', 'build', 'src')
bustools = op.join(config['working_dir'], 'software', 'bustools', 'build', 'src')
shell.prefix('export PATH=' + kallisto + ':' + bustools + ":$PATH;")

try:
    os.makedirs(op.join(config['working_dir'], 'logs'), exist_ok=True)
    os.makedirs(op.join(config['working_dir'], 'benchmarks'), exist_ok=True)
except OSError:
    pass

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

rule all:
    input:
        expand(op.join(config['working_dir'], '{aligner}',  '{sample}', 'descriptive_report.html'),
               aligner = get_aligners(),
               sample = get_sample_names()),
        _sampletag_reports,
        _extra_targets
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

rule compile_kallisto:
    output:
        op.join(config['working_dir'], 'software', 'kallisto', 'build', 'src', 'kallisto')
    log:
        op.join(config['working_dir'], 'logs', 'kallisto_install.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'kallisto_install.txt')
    params:
        soft = op.join(config['working_dir'], 'software')
    threads:
        1 # min(5, workflow.cores)
    shell:
        """
        mkdir -p {params.soft}
        cd {params.soft}
        
        rm -rf kallisto
        # kallisto
        git clone https://github.com/pachterlab/kallisto.git --depth 1
        cd kallisto 
        git log | head &> {log}
        mkdir build
        cd build
        cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME &>> {log}
        make -j {threads} &>> {log}
        """

rule compile_bustools:
    output:
        op.join(config['working_dir'], 'software', 'bustools', 'build', 'src', 'bustools')
    log:
        op.join(config['working_dir'], 'logs', 'bustools_install.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'bustools_install.txt')
    params:
        soft = op.join(config['working_dir'], 'software')
    threads:
        min(5, workflow.cores)
    shell:
        """
        mkdir -p {params.soft}
        cd {params.soft}
        rm -rf bustools
        
        # bustools
       
        git clone https://github.com/BUStools/bustools.git --depth 1
        cd bustools
        git log | head &> {log}
        mkdir build
        cd build
        cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME &>> {log}
        make -j {threads} &>> {log}
        """

        
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
        star = config['STAR'],
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

    ({params.star} --runThreadN {params.nthreads} \
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
                                  'filtered', 'matrix.mtx')
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
        STAR = config['STAR'],
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

   {params.STAR} --runThreadN {params.threads} \
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

## yes the log is considered an output - to pass as a flag
## R_LIBS are conda's if run in conda, but /home/rock/R_LIBs if run in docker, and user's if run directly
rule install_r_deps:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        script = op.join(config['repo_path'], 'src', 'installs.R')
    output:
        log = op.join(config['working_dir'], 'logs', 'installs.log')
    params:
        working_dir = config['working_dir'],
        Rbin = config['Rbin'],
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_install.txt')
    shell:
        """
        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} &> {output.log}
         """

# todo fixme so it gets the filtered mtx
rule generate_sce_starsolo:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        filtered  = op.join(config['working_dir'], 'starsolo', '{sample}', 'Solo.out',
                               'Gene', 'filtered', 'matrix.mtx'),
        bam = op.join(config['working_dir'], 'starsolo', '{sample}', 'Aligned.sortedByCoord.out.bam'),
        script = op.join(config['repo_path'], 'src', 'generate_sce_star.R'),
        installs = op.join(config['working_dir'], 'logs', 'installs.log')
    output:
        sce = op.join(config['working_dir'], 'starsolo', '{sample}', '{sample}_starsolo_sce.rds')
    params:
        align_path = op.join(config['working_dir']),
        working_dir = config['working_dir'],
        sample = lambda wildcards: wildcards.sample,
        Rbin = config['Rbin']
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_star.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_star.txt')
    shell:
        """
        ## this is unrelated to the bamgeneration; fixes starsolo's default permissions
        chmod -R ug+rwX $(dirname {input.bam})

        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --output_fn {output.sce} &> {log}
        """

# todo fixme so it gets the filtered mtx
rule generate_sce_kallisto:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        flag = op.join(config['working_dir'], 'bustools', '{sample}', 'output.mtx'),
        # gtf = config['gtf'],
        script = op.join(config['repo_path'], 'src', 'generate_sce_kallisto.R'),
        installs = op.join(config['working_dir'], 'logs', 'installs.log')
    output:
        sce = op.join(config['working_dir'], 'kallisto', '{sample}', '{sample}_kallisto_sce.rds')
    params:
        working_dir = config['working_dir'],
        sample = lambda wildcards: wildcards.sample,
        Rbin = config['Rbin']
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_kallisto.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_kallisto.txt')
    shell:
        """
        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
             --output_fn {output.sce} &> {log}
        """

rule generate_sce_alevin:
    conda:
        op.join('envs', 'all_in_one.yaml')
    input:
        flag = op.join(config['working_dir'], 'alevin', '{sample}', 'alevin', 'quants_mat.gz'),
        script = op.join(config['repo_path'], 'src', 'generate_sce_alevin.R'),
        installs = op.join(config['working_dir'], 'logs', 'installs.log')
    output:
        sce = op.join(config['working_dir'], 'alevin', '{sample}', '{sample}_alevin_sce.rds')
    params:
        working_dir = config['working_dir'],
        sample = lambda wildcards: wildcards.sample,
        Rbin = config['Rbin']
    log:
        op.join(config['working_dir'], 'logs', 'r_sce_generation_{sample}_alevin.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', 'r_sce_generation_{sample}_alevin.txt')
    shell:
        """
        {params.Rbin} -q --no-save --no-restore --slave \
             -f {input.script} --args \
             --sample {wildcards.sample} \
             --working_dir {params.working_dir} \
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
#         Rbin = config['Rbin']
#     shell:
#         """
#         {params.Rbin} -q --no-save --no-restore --slave \
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
        installs = op.join(config['working_dir'], 'logs', 'installs.log'),
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
        Rbin = config['Rbin']
    shell:
        """
        cd {params.working_dir}
        mkdir -p {params.path}
        {params.Rbin} --vanilla -e 'rmarkdown::render(\"{input.script}\", 
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

## conda recipe is broken        
rule kallisto_index:
    # conda:
    #     op.join('envs', 'kallisto.yaml')
    input:
        transcriptome = config['transcriptome'],
        kal = op.join(config['working_dir'], 'software', 'kallisto', 'build', 'src', 'kallisto'),
        bus = op.join(config['working_dir'], 'software', 'bustools', 'build', 'src', 'bustools')  
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

## conda recipe is broken 
rule kallisto_bus:
    # conda:
    #     op.join('envs', 'kallisto.yaml')
    input:
        transcriptome = config['transcriptome'],
        # transcriptome = op.join(config['working_dir'], 'data', 'index', 'salmon', 'transcriptome.fa'),        
        # cdna = lambda wildcards: get_cdna_by_name(wildcards.sample),
        standardized_cdna = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cdna.fq.gz"),
        standardized_cb_umi = op.join(config['working_dir'], 'data', 'fastq', "{sample}_standardized_cb_umi.fq.gz"),
        kallisto_index = op.join(config['working_dir'], 'data', 'index', 'kallisto', 'kallisto.index'),
        kal = op.join(config['working_dir'], 'software', 'kallisto', 'build', 'src', 'kallisto'),
        bus = op.join(config['working_dir'], 'software', 'bustools', 'build', 'src', 'bustools')
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
    input:
        bus    = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.bus'),
        btools = op.join(config['working_dir'], 'software', 'bustools', 'build', 'src', 'bustools')
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

rule bustools_count:
    # conda:
    #     op.join('envs', 'kallisto.yaml')
    input:
        txp2gene   = op.join(config['working_dir'], 'data', 'index', 'salmon', 'txp2gene'),
        matrix_ec  = op.join(config['working_dir'], 'kallisto', '{sample}', 'matrix.ec'),
        transcripts = op.join(config['working_dir'], 'kallisto', '{sample}', 'transcripts.txt'),
        bus        = op.join(config['working_dir'], 'kallisto', '{sample}', 'output.sorted.bus'),
        kal        = op.join(config['working_dir'], 'software', 'kallisto', 'build', 'src', 'kallisto'),
        btools     = op.join(config['working_dir'], 'software', 'bustools', 'build', 'src', 'bustools')
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
        installs = op.join(config['working_dir'], 'logs', 'installs.log')
    output:
        html = op.join(config['working_dir'], 'sampletags', '{sample}', 'sampletag_report.html')
    log:
        op.join(config['working_dir'], 'logs', '{sample}_sampletag_report.log')
    benchmark:
        op.join(config['working_dir'], 'benchmarks', '{sample}_sampletag_report.txt')
    params:
        Rbin = config['Rbin'],
        assignments = lambda wildcards: (
            op.join(config['working_dir'], 'simulate', 'sampletag_assignments.txt')
            if config.get('use_simulated', False) else ''
        )
    shell:
        """
        {params.Rbin} --vanilla -e \
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
            zcat {input.transcriptome} | awk '/^>/ {{sub(/\..*/, "", $1)}} {{print}}' > {output.deversioned_fasta}

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
        op.join(config['working_dir'], 'benchmarks', 'alevin_{sample}_align.log')
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
#         STAR = config['STAR'],
#         # num_cells = get_expected_cells_by_name("{sample}"),
#         tmp = op.join(config['working_dir'], 'tmp_tasseq_{sample}'),
#         maxmem = config['max_mem_mb'] * 1024 * 1024,
#         sjdbOverhang = config['sjdbOverhang'],
#         soloCellFilter = config['soloCellFilter']
#     shell:
#         """
#         rm -rf {params.tmp}
#         mkdir -p {params.path} 

#         {params.STAR} --runThreadN {params.threads} \
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
