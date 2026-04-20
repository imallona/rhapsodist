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
def _raw_cbumi_path(name):
    for s in config['samples']:
        if s['name'] == name:
            if 'cb_umi_fq' in s['uses']:
                return s['uses']['cb_umi_fq']
            elif s['uses'].get('sra_run'):
                return op.join(config['working_dir'], 'data', 'fastq', 'sra', name, name + '_R1.fastq.gz')

def _raw_cdna_path(name):
    for s in config['samples']:
        if s['name'] == name:
            if 'cdna_fq' in s['uses']:
                return s['uses']['cdna_fq']
            elif s['uses'].get('sra_run'):
                return op.join(config['working_dir'], 'data', 'fastq', 'sra', name, name + '_R2.fastq.gz')

def get_downsample_fraction(name):
    """Return the downsample fraction in [0, 1]. Per-sample 'downsample' overrides
    the global 'downsample' config key. Value is a percentage (0.01 to 100);
    default is 100 (keep all reads). Values <= 0 or > 100 raise."""
    uses = _sample_uses(name) or {}
    pct = uses.get('downsample', config.get('downsample', 100))
    pct = float(pct)
    if pct <= 0 or pct > 100:
        raise ValueError(
            f"Sample '{name}': downsample must be in (0, 100], got {pct}"
        )
    return pct / 100.0

def _downsampled_fastq(name, mate):
    return op.join(config['working_dir'], 'data', 'fastq', 'downsampled',
                   f'{name}_{mate}.fastq.gz')

def get_cbumi_by_name(name):
    if get_downsample_fraction(name) < 1.0:
        return _downsampled_fastq(name, 'R1')
    return _raw_cbumi_path(name)

def get_cdna_by_name(name):
    if get_downsample_fraction(name) < 1.0:
        return _downsampled_fastq(name, 'R2')
    return _raw_cdna_path(name)

def get_expected_cells_by_name(name):
    for i in range(len(config['samples'])):
        if config['samples'][i]['name'] == name:
             return(config['samples'][i]['uses']['expected_cells'])

def detect_bead_version(cb_umi_path, n_reads=10000):
    """Detect bead class from the first n_reads sequence reads of an R1 fastq.

    Returns 'v1', 'enhanced', or 'unknown'.
    v1 beads have linker ACTGGCCTGCGA at read positions 10-21 (0-indexed 9:21).
    Enhanced beads have 4bp linker GTGA at positions 10-13, shifted by 0-3 bp
    stagger (so the motif can sit anywhere in positions 10-16).
    """
    import gzip as _gzip
    count = 0
    v1_count = 0
    enh_count = 0
    opener = _gzip.open if str(cb_umi_path).endswith('.gz') else open
    with opener(cb_umi_path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            seq = line.strip()
            count += 1
            if count > n_reads:
                break
            if len(seq) >= 21 and seq[9:21] == 'ACTGGCCTGCGA':
                v1_count += 1
            else:
                for offset in range(4):
                    if len(seq) >= 13 + offset and seq[9 + offset:13 + offset] == 'GTGA':
                        enh_count += 1
                        break
    if count == 0:
        return 'unknown'
    v1_frac = v1_count / count
    enh_frac = enh_count / count
    if v1_frac >= enh_frac and v1_frac > 0.1:
        return 'v1'
    if enh_frac > 0.1:
        return 'enhanced'
    return 'unknown'

def _sample_uses(name):
    for s in config['samples']:
        if s['name'] == name:
            return s['uses']
    return None

def get_declared_allowedlist(name):
    """Return declared allowedlist ('96' or '384') or None. QC/printing only."""
    uses = _sample_uses(name) or {}
    v = uses.get('allowedlist')
    if v is None:
        return None
    v = str(v)
    if v not in ('96', '384'):
        raise ValueError(
            f"Sample '{name}': allowedlist must be 96 or 384, got '{v}'"
        )
    return v

def _yesno(name, field, v):
    """Coerce a config value to 'yes'/'no'. Accepts bool (True/False) from YAML's
    unquoted yes/no, or the literal strings 'yes'/'no' (case-insensitive)."""
    if isinstance(v, bool):
        return 'yes' if v else 'no'
    s = str(v).lower()
    if s in ('yes', 'true'):
        return 'yes'
    if s in ('no', 'false'):
        return 'no'
    raise ValueError(
        f"Sample '{name}': {field} must be yes or no, got '{v}'"
    )

def get_declared_diversity_insets(name):
    """Return declared diversity_insets ('yes' or 'no') or None. QC/printing only."""
    uses = _sample_uses(name) or {}
    v = uses.get('diversity_insets')
    if v is None:
        return None
    return _yesno(name, 'diversity_insets', v)

def get_bead_type_by_name(name):
    """Return internal bead class ('v1', 'enhanced', 'enhanced_v2') derived from
    declared allowedlist + diversity_insets. Defaults to 'enhanced_v2' when both
    are absent. The actual chemistry is still checked by detect_bead_version and
    reported in the validation log."""
    al = get_declared_allowedlist(name)
    di = get_declared_diversity_insets(name)
    if al is None and di is None:
        return 'enhanced_v2'
    al = al or '384'
    di = di or 'yes'
    if al == '96' and di == 'no':
        return 'v1'
    if al == '96' and di == 'yes':
        return 'enhanced'
    if al == '384' and di == 'yes':
        return 'enhanced_v2'
    raise ValueError(
        f"Sample '{name}': unsupported combination allowedlist={al}, "
        f"diversity_insets={di} (valid: 96/no, 96/yes, 384/yes)"
    )

def get_barcode_whitelist_by_name(name):
    """Return whitelist directory name ('96x3' or '384x3') from declared allowedlist."""
    al = get_declared_allowedlist(name)
    if al is None:
        return '384x3' if get_bead_type_by_name(name) == 'enhanced_v2' else '96x3'
    return al + 'x3'

def get_guide_url_by_name(name):
    """Return NCBI FTP URL for the guide assignment CSV, or None if not configured."""
    for s in config['samples']:
        if s['name'] == name:
            gsm = s['uses'].get('guide_gsm')
            wta = s['uses'].get('guide_wta')
            if gsm and wta:
                prefix = gsm[:-3] + 'nnn'
                fn = f"{gsm}_guides_dialout_{wta}_umi_counts_anno.csv.gz"
                return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix}/{gsm}/suppl/{fn}"
    return None

def get_use_sampletags(name):
    """True when sampletag demux runs for this sample. Per-sample 'use_sampletags'
    (yes/no) takes precedence; otherwise follows the global 'skip_sampletags' flag."""
    uses = _sample_uses(name) or {}
    v = uses.get('use_sampletags')
    if v is not None:
        return _yesno(name, 'use_sampletags', v) == 'yes'
    return not config.get('skip_sampletags', False)

def get_species_by_name(name):
    """Return declared species ('human' or 'mouse') or None when not declared.
    When use_sampletags is yes, species is mandatory and this raises if absent."""
    uses = _sample_uses(name) or {}
    species = uses.get('species')
    if species is None:
        if get_use_sampletags(name):
            raise ValueError(
                f"Sample '{name}': species is required when use_sampletags is yes"
            )
        return None
    if species not in ('mouse', 'human'):
        raise ValueError(
            f"Sample '{name}': unknown species '{species}' (expected mouse or human)"
        )
    return species

def samples_with_sampletags():
    return [s for s in get_sample_names() if get_use_sampletags(s)]

def validate_sampletag_config():
    """Raise early if any sample has use_sampletags=yes without species."""
    for name in get_sample_names():
        if get_use_sampletags(name):
            get_species_by_name(name)
            
             
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
    """Derive Sample_Tags_Version from per-sample override, global config, or species.
    Returns None when nothing is declared (sample runs without sampletag demux)."""
    v = _sbg_uses(name, 'sbg_sample_tags_version')
    if v:
        return v
    return get_species_by_name(name)

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
