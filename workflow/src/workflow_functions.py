#!/usr/bin/env python
##
## Functions process rock/roi data (general method)
##
## Started 11th Oct 2023
##
## Izaskun Mallona
## GPLv3

def get_sample_names():
    return([x['name'] for x in config['samples']])

def get_aligners():
    return(config['aligner'])

## name means sample name, everywhere
def get_cbumi_by_name(name):
    for s in config['samples']:
        if s['name'] == name:
            if 'cb_umi_fq' in s['uses']:
                return s['uses']['cb_umi_fq']
            elif 'sra_run' in s['uses']:
                return op.join(config['working_dir'], 'data', 'fastq', 'sra', name, name + '_R1.fastq.gz')

def get_cdna_by_name(name):
    for s in config['samples']:
        if s['name'] == name:
            if 'cdna_fq' in s['uses']:
                return s['uses']['cdna_fq']
            elif 'sra_run' in s['uses']:
                return op.join(config['working_dir'], 'data', 'fastq', 'sra', name, name + '_R2.fastq.gz')

def get_expected_cells_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             return(config['samples'][i]['uses']['expected_cells'])

def get_bead_type_by_name(name):
    """Return bead_version for a sample: 'v1', 'enhanced', or 'enhanced_v2' (default)."""
    for s in config['samples']:
        if s['name'] == name:
            return s['uses'].get('bead_version', 'enhanced_v2')
    return 'enhanced_v2'

def get_barcode_whitelist_by_name(name):
    """Derive whitelist directory name from bead_version."""
    bv = get_bead_type_by_name(name)
    if bv == 'enhanced_v2':
        return '384x3'
    return '96x3'

def get_species_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             species = config['samples'][i]['uses']['species']
             if species in ['mouse', 'human']:
                 return(species)
             else:
                 raise ValueError('Unknown species (not mouse nor human), it was reported ' + species)
             
# def get_sampletags_fasta_by_name(name):
#     for i in range(len(config['samples'])):
#         if config['samples'][i]['name'] == name:
#             species = config['samples'][i]['uses']['species']
#             return op.join('data', 'sampletags', species + '_sampletags.fa')
            
             
def get_chromosomes(wildcards):
    # with open(op.join(config['working_dir'], 'data', 'chrom.sizes')) as fh:
    # with open(chromsizes_fn) as fh:
    fn = checkpoints.retrieve_genome_sizes.get(**wildcards).output[0]
    with open(fn) as fh:
        return(list(line.strip().split('\t')[0] for line in fh))

# def list_by_chr_dedup_bams(wildcards):
#     chroms = get_chromosomes(wildcards)
#     return(chrom + '_cb_umi_deduped.bam' for chrom in chroms)

## ── SBG / BD Rhapsody official pipeline helpers ──────────────────────────────
## Each function checks the per-sample uses: block first, then falls back to
## the top-level config.  This lets users specify SBG options once globally or
## override them per sample.

def _sbg_uses(name, key):
    """Return per-sample SBG config key, or global fallback, or None."""
    for s in config['samples']:
        if s['name'] == name:
            return s['uses'].get(key) or config.get(key)
    return config.get(key)

def get_sbg_cwl_by_name(name):
    return _sbg_uses(name, 'sbg_cwl')

def get_sbg_reference_by_name(name):
    return _sbg_uses(name, 'sbg_reference_archive')

def get_sbg_sample_tags_version_by_name(name):
    """Derive Sample_Tags_Version: per-sample override → global config → species."""
    v = _sbg_uses(name, 'sbg_sample_tags_version')
    if v:
        return v
    return get_species_by_name(name)   # 'human' or 'mouse'

def get_sbg_bead_version_by_name(name):
    explicit = _sbg_uses(name, 'sbg_bead_version') or config.get('sbg_bead_version')
    if explicit:
        return explicit
    bv = get_bead_type_by_name(name)
    mapping = {'v1': 'Multiplex', 'enhanced': 'Enh', 'enhanced_v2': 'EnhV2'}
    return mapping.get(bv, 'Enh')

def get_sbg_reference_url_by_name(name):
    return _sbg_uses(name, 'sbg_reference_url') or config.get('sbg_reference_url')

## bd offers a couple of sets of whitelists, so we fetch the right one according to the config.yaml file
def symlink_whitelist(sample):
    os.makedirs(op.join(config['working_dir'], 'starsolo'), exist_ok = True)
                
    if get_barcode_whitelist_by_name(name = sample) == '96x3':
        for x in ['BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt']:
            try:
                os.symlink(src = op.join(workflow.basedir, 'data', 'whitelist_96x3', x),
                           dst = op.join(config['working_dir'], 'starsolo', sample, 'whitelists', x))
            except FileExistsError:
                continue
    elif get_barcode_whitelist_by_name(name = sample) == '384x3':
        for x in ['BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt']:
            try:
                os.symlink(src = op.join(workflow.basedir, 'data', 'whitelist_384x3', x),
                           dst = op.join(config['working_dir'], 'starsolo', sample, 'whitelists', x))
            except FileExistsError:
                continue
