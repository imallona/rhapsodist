import gzip
import os
import os.path as op

import pytest

# workflow_functions.py is designed to be exec'd by Snakemake (include:), so it
# has no imports of its own. Import it as a normal module and inject the globals
# that Snakemake would normally provide.
import workflow.src.workflow_functions as wf

wf.os = os
wf.op = op


SAMPLE_A = {
    'name': 'sampleA',
    'uses': {
        'cb_umi_fq': '/data/sampleA_R1.fq.gz',
        'cdna_fq': '/data/sampleA_R2.fq.gz',
        'allowedlist': 96,
        'diversity_insets': 'yes',
        'use_sampletags': 'yes',
        'species': 'human',
    },
}

SAMPLE_B = {
    'name': 'sampleB',
    'uses': {
        'cb_umi_fq': '/data/sampleB_R1.fq.gz',
        'cdna_fq': '/data/sampleB_R2.fq.gz',
        'allowedlist': 384,
        'diversity_insets': 'yes',
        'use_sampletags': 'yes',
        'species': 'mouse',
        'sbg_cwl': '/per_sample/pipeline.cwl',
    },
}


def make_config(**extra):
    cfg = {
        'aligner': ['starsolo', 'kallisto'],
        'samples': [
            {'name': s['name'], 'uses': dict(s['uses'])}
            for s in (SAMPLE_A, SAMPLE_B)
        ],
    }
    cfg.update(extra)
    return cfg


@pytest.fixture(autouse=True)
def reset_config():
    wf.config = make_config()
    yield
    wf.config = {}


def test_get_sample_names():
    assert wf.get_sample_names() == ['sampleA', 'sampleB']


def test_get_aligners():
    assert wf.get_aligners() == ['starsolo', 'kallisto']


def test_get_cbumi_by_name():
    assert wf.get_cbumi_by_name('sampleA') == '/data/sampleA_R1.fq.gz'


def test_get_cdna_by_name():
    assert wf.get_cdna_by_name('sampleA') == '/data/sampleA_R2.fq.gz'


def test_get_cbumi_by_name_sra(tmp_path):
    wf.config['working_dir'] = str(tmp_path)
    wf.config['samples'][0]['uses'].pop('cb_umi_fq')
    wf.config['samples'][0]['uses']['sra_run'] = 'SRR123'
    result = wf.get_cbumi_by_name('sampleA')
    assert result.endswith('sampleA_R1.fastq.gz')


def test_get_cdna_by_name_sra(tmp_path):
    wf.config['working_dir'] = str(tmp_path)
    wf.config['samples'][0]['uses'].pop('cdna_fq')
    wf.config['samples'][0]['uses']['sra_run'] = 'SRR123'
    result = wf.get_cdna_by_name('sampleA')
    assert result.endswith('sampleA_R2.fastq.gz')


def test_get_expected_cells_by_name():
    wf.config['samples'][0]['uses']['expected_cells'] = 5000
    assert wf.get_expected_cells_by_name('sampleA') == 5000


## allowedlist / diversity_insets / bead type derivation -----------------------

def test_get_declared_allowedlist_96():
    assert wf.get_declared_allowedlist('sampleA') == '96'


def test_get_declared_allowedlist_384():
    assert wf.get_declared_allowedlist('sampleB') == '384'


def test_get_declared_allowedlist_absent():
    wf.config['samples'][0]['uses'].pop('allowedlist')
    assert wf.get_declared_allowedlist('sampleA') is None


def test_get_declared_allowedlist_invalid():
    wf.config['samples'][0]['uses']['allowedlist'] = 42
    with pytest.raises(ValueError):
        wf.get_declared_allowedlist('sampleA')


def test_get_declared_diversity_insets_yes():
    assert wf.get_declared_diversity_insets('sampleA') == 'yes'


def test_get_declared_diversity_insets_no():
    wf.config['samples'][0]['uses']['diversity_insets'] = 'no'
    assert wf.get_declared_diversity_insets('sampleA') == 'no'


def test_get_declared_diversity_insets_absent():
    wf.config['samples'][0]['uses'].pop('diversity_insets')
    assert wf.get_declared_diversity_insets('sampleA') is None


def test_get_declared_diversity_insets_invalid():
    wf.config['samples'][0]['uses']['diversity_insets'] = 'maybe'
    with pytest.raises(ValueError):
        wf.get_declared_diversity_insets('sampleA')


def test_get_bead_type_96_no_is_v1():
    wf.config['samples'][0]['uses']['allowedlist'] = 96
    wf.config['samples'][0]['uses']['diversity_insets'] = 'no'
    assert wf.get_bead_type_by_name('sampleA') == 'v1'


def test_get_bead_type_96_yes_is_enhanced():
    wf.config['samples'][0]['uses']['allowedlist'] = 96
    wf.config['samples'][0]['uses']['diversity_insets'] = 'yes'
    assert wf.get_bead_type_by_name('sampleA') == 'enhanced'


def test_get_bead_type_384_yes_is_enhanced_v2():
    wf.config['samples'][0]['uses']['allowedlist'] = 384
    wf.config['samples'][0]['uses']['diversity_insets'] = 'yes'
    assert wf.get_bead_type_by_name('sampleA') == 'enhanced_v2'


def test_get_bead_type_384_no_is_invalid():
    wf.config['samples'][0]['uses']['allowedlist'] = 384
    wf.config['samples'][0]['uses']['diversity_insets'] = 'no'
    with pytest.raises(ValueError):
        wf.get_bead_type_by_name('sampleA')


def test_get_bead_type_both_absent_defaults_enhanced_v2():
    wf.config['samples'][0]['uses'].pop('allowedlist')
    wf.config['samples'][0]['uses'].pop('diversity_insets')
    assert wf.get_bead_type_by_name('sampleA') == 'enhanced_v2'


def test_get_bead_type_fallback_unknown_sample():
    wf.config['samples'] = []
    assert wf.get_bead_type_by_name('nonexistent') == 'enhanced_v2'


def test_get_barcode_whitelist_96x3():
    assert wf.get_barcode_whitelist_by_name('sampleA') == '96x3'


def test_get_barcode_whitelist_384x3():
    assert wf.get_barcode_whitelist_by_name('sampleB') == '384x3'


def test_get_barcode_whitelist_defaults_to_384x3_when_absent():
    wf.config['samples'][0]['uses'].pop('allowedlist')
    wf.config['samples'][0]['uses'].pop('diversity_insets')
    assert wf.get_barcode_whitelist_by_name('sampleA') == '384x3'


## use_sampletags and species -------------------------------------------------

def test_get_use_sampletags_yes():
    assert wf.get_use_sampletags('sampleA') is True


def test_get_use_sampletags_no_per_sample():
    wf.config['samples'][0]['uses']['use_sampletags'] = 'no'
    wf.config['samples'][0]['uses'].pop('species')
    assert wf.get_use_sampletags('sampleA') is False


def test_get_use_sampletags_follows_global_skip_when_absent():
    wf.config['samples'][0]['uses'].pop('use_sampletags')
    wf.config['skip_sampletags'] = True
    assert wf.get_use_sampletags('sampleA') is False


def test_get_use_sampletags_follows_global_not_skip_when_absent():
    wf.config['samples'][0]['uses'].pop('use_sampletags')
    wf.config['skip_sampletags'] = False
    assert wf.get_use_sampletags('sampleA') is True


def test_get_use_sampletags_invalid():
    wf.config['samples'][0]['uses']['use_sampletags'] = 'maybe'
    with pytest.raises(ValueError):
        wf.get_use_sampletags('sampleA')


def test_get_species_human():
    assert wf.get_species_by_name('sampleA') == 'human'


def test_get_species_mouse():
    assert wf.get_species_by_name('sampleB') == 'mouse'


def test_get_species_absent_when_use_sampletags_no():
    wf.config['samples'][0]['uses']['use_sampletags'] = 'no'
    wf.config['samples'][0]['uses'].pop('species')
    assert wf.get_species_by_name('sampleA') is None


def test_get_species_required_when_use_sampletags_yes():
    wf.config['samples'][0]['uses']['use_sampletags'] = 'yes'
    wf.config['samples'][0]['uses'].pop('species')
    with pytest.raises(ValueError):
        wf.get_species_by_name('sampleA')


def test_get_species_invalid_raises():
    wf.config['samples'][0]['uses']['species'] = 'zebrafish'
    with pytest.raises(ValueError):
        wf.get_species_by_name('sampleA')


def test_samples_with_sampletags_all():
    assert wf.samples_with_sampletags() == ['sampleA', 'sampleB']


def test_samples_with_sampletags_none():
    for s in wf.config['samples']:
        s['uses']['use_sampletags'] = 'no'
        s['uses'].pop('species', None)
    assert wf.samples_with_sampletags() == []


def test_samples_with_sampletags_mixed():
    wf.config['samples'][1]['uses']['use_sampletags'] = 'no'
    wf.config['samples'][1]['uses'].pop('species')
    assert wf.samples_with_sampletags() == ['sampleA']


def test_validate_sampletag_config_passes():
    wf.validate_sampletag_config()  # should not raise


def test_validate_sampletag_config_catches_missing_species():
    wf.config['samples'][0]['uses']['use_sampletags'] = 'yes'
    wf.config['samples'][0]['uses'].pop('species')
    with pytest.raises(ValueError):
        wf.validate_sampletag_config()


## downsample -----------------------------------------------------------------

def test_get_downsample_fraction_default_is_one():
    assert wf.get_downsample_fraction('sampleA') == 1.0


def test_get_downsample_fraction_global():
    wf.config['downsample'] = 50
    assert wf.get_downsample_fraction('sampleA') == 0.5


def test_get_downsample_fraction_per_sample_overrides_global():
    wf.config['downsample'] = 50
    wf.config['samples'][0]['uses']['downsample'] = 10
    assert wf.get_downsample_fraction('sampleA') == 0.1


def test_get_downsample_fraction_accepts_fractional_percentage():
    wf.config['samples'][0]['uses']['downsample'] = 0.01
    assert wf.get_downsample_fraction('sampleA') == pytest.approx(0.0001)


def test_get_downsample_fraction_rejects_zero():
    wf.config['samples'][0]['uses']['downsample'] = 0
    with pytest.raises(ValueError):
        wf.get_downsample_fraction('sampleA')


def test_get_downsample_fraction_rejects_over_100():
    wf.config['samples'][0]['uses']['downsample'] = 150
    with pytest.raises(ValueError):
        wf.get_downsample_fraction('sampleA')


def test_get_cbumi_uses_downsampled_path_when_downsample_enabled(tmp_path):
    wf.config['working_dir'] = str(tmp_path)
    wf.config['samples'][0]['uses']['downsample'] = 10
    result = wf.get_cbumi_by_name('sampleA')
    assert 'downsampled' in result
    assert result.endswith('sampleA_R1.fastq.gz')


def test_get_cdna_uses_downsampled_path_when_downsample_enabled(tmp_path):
    wf.config['working_dir'] = str(tmp_path)
    wf.config['samples'][0]['uses']['downsample'] = 10
    result = wf.get_cdna_by_name('sampleA')
    assert 'downsampled' in result
    assert result.endswith('sampleA_R2.fastq.gz')


def test_get_cbumi_uses_raw_path_when_downsample_full():
    assert wf.get_cbumi_by_name('sampleA') == '/data/sampleA_R1.fq.gz'


## multiple input fastqs per sample ------------------------------------------

def test_single_fastq_as_list_returns_the_path():
    wf.config['samples'][0]['uses']['cb_umi_fq'] = ['/data/sampleA_R1.fq.gz']
    wf.config['samples'][0]['uses']['cdna_fq'] = ['/data/sampleA_R2.fq.gz']
    assert wf.get_cbumi_by_name('sampleA') == '/data/sampleA_R1.fq.gz'
    assert wf.get_cdna_by_name('sampleA') == '/data/sampleA_R2.fq.gz'


def test_multiple_fastqs_return_combined_path(tmp_path):
    wf.config['working_dir'] = str(tmp_path)
    wf.config['samples'][0]['uses']['cb_umi_fq'] = ['/data/a_R1.fq.gz', '/data/b_R1.fq.gz']
    wf.config['samples'][0]['uses']['cdna_fq'] = ['/data/a_R2.fq.gz', '/data/b_R2.fq.gz']
    r1 = wf.get_cbumi_by_name('sampleA')
    r2 = wf.get_cdna_by_name('sampleA')
    assert 'combined' in r1 and r1.endswith('sampleA_R1.fastq.gz')
    assert 'combined' in r2 and r2.endswith('sampleA_R2.fastq.gz')


def test_get_fastq_inputs_return_lists():
    wf.config['samples'][0]['uses']['cb_umi_fq'] = ['/data/a_R1.fq.gz', '/data/b_R1.fq.gz']
    assert wf.get_cbumi_inputs('sampleA') == ['/data/a_R1.fq.gz', '/data/b_R1.fq.gz']
    assert wf.get_cdna_inputs('sampleA') == ['/data/sampleA_R2.fq.gz']


def test_get_fastq_inputs_none_for_sra():
    wf.config['samples'][0]['uses'].pop('cb_umi_fq')
    assert wf.get_cbumi_inputs('sampleA') is None


def test_validate_fastq_lists_passes_on_matching_lengths():
    wf.config['samples'][0]['uses']['cb_umi_fq'] = ['/data/a_R1.fq.gz', '/data/b_R1.fq.gz']
    wf.config['samples'][0]['uses']['cdna_fq'] = ['/data/a_R2.fq.gz', '/data/b_R2.fq.gz']
    wf.validate_fastq_lists()  # should not raise


def test_validate_fastq_lists_catches_length_mismatch():
    wf.config['samples'][0]['uses']['cb_umi_fq'] = ['/data/a_R1.fq.gz', '/data/b_R1.fq.gz']
    wf.config['samples'][0]['uses']['cdna_fq'] = ['/data/a_R2.fq.gz']
    with pytest.raises(ValueError):
        wf.validate_fastq_lists()


## guide URL, SBG, etc. -------------------------------------------------------

def test_get_guide_url_returns_none_when_no_gsm():
    assert wf.get_guide_url_by_name('sampleA') is None


def test_get_guide_url_constructs_url():
    wf.config['samples'][0]['uses']['guide_gsm'] = 'GSM7500353'
    wf.config['samples'][0]['uses']['guide_wta'] = 'WTA16'
    url = wf.get_guide_url_by_name('sampleA')
    assert 'GSM7500353' in url
    assert 'WTA16' in url
    assert url.startswith('https://')


def test_get_sbg_cwl_by_name():
    wf.config['samples'][0]['uses']['sbg_cwl'] = '/path/to/pipeline.cwl'
    assert wf.get_sbg_cwl_by_name('sampleA') == '/path/to/pipeline.cwl'


def test_sbg_uses_per_sample_key():
    assert wf._sbg_uses('sampleB', 'sbg_cwl') == '/per_sample/pipeline.cwl'


def test_sbg_uses_global_fallback():
    wf.config['sbg_cwl'] = '/global/pipeline.cwl'
    assert wf._sbg_uses('sampleA', 'sbg_cwl') == '/global/pipeline.cwl'


def test_sbg_uses_returns_none_when_absent():
    assert wf._sbg_uses('sampleA', 'sbg_cwl') is None


def test_get_sbg_bead_version_384x3_gives_enhv2():
    assert wf.get_sbg_bead_version_by_name('sampleB') == 'EnhV2'


def test_get_sbg_bead_version_96x3_gives_enh():
    assert wf.get_sbg_bead_version_by_name('sampleA') == 'Enh'


def test_get_sbg_bead_version_explicit_override():
    wf.config['samples'][0]['uses']['sbg_bead_version'] = 'v1'
    assert wf.get_sbg_bead_version_by_name('sampleA') == 'v1'


def test_get_sbg_sample_tags_version_falls_back_to_species():
    assert wf.get_sbg_sample_tags_version_by_name('sampleA') == 'human'


def test_get_sbg_sample_tags_version_explicit():
    wf.config['samples'][0]['uses']['sbg_sample_tags_version'] = 'mouse'
    assert wf.get_sbg_sample_tags_version_by_name('sampleA') == 'mouse'


def test_get_sbg_sample_tags_version_none_when_no_species():
    wf.config['samples'][0]['uses']['use_sampletags'] = 'no'
    wf.config['samples'][0]['uses'].pop('species')
    assert wf.get_sbg_sample_tags_version_by_name('sampleA') is None


def test_get_sbg_reference_url_none_when_absent():
    assert wf.get_sbg_reference_url_by_name('sampleA') is None


def test_get_sbg_reference_url_from_global():
    wf.config['sbg_reference_url'] = 'http://example.com/ref.tar.gz'
    assert wf.get_sbg_reference_url_by_name('sampleA') == 'http://example.com/ref.tar.gz'


def test_get_sbg_reference_none_when_absent():
    assert wf.get_sbg_reference_by_name('sampleA') is None


def test_get_sbg_reference_from_global():
    wf.config['sbg_reference_archive'] = '/data/ref.tar.gz'
    assert wf.get_sbg_reference_by_name('sampleA') == '/data/ref.tar.gz'


## bead class detection -------------------------------------------------------

def _write_fastq(path, sequences):
    with open(path, 'w') as fh:
        for i, seq in enumerate(sequences):
            fh.write(f'@read{i}\n{seq}\n+\n{"F" * len(seq)}\n')


def _write_fastq_gz(path, sequences):
    with gzip.open(path, 'wt') as fh:
        for i, seq in enumerate(sequences):
            fh.write(f'@read{i}\n{seq}\n+\n{"F" * len(seq)}\n')


V1_READ = 'A' * 9 + 'ACTGGCCTGCGA' + 'C' * 9 + 'GGTAGCGGTGACA' + 'G' * 9 + 'T' * 8
ENH_READ = 'A' * 9 + 'GTGA' + 'C' * 9 + 'GACA' + 'G' * 9 + 'T' * 8
ENH_READ_VB1 = 'N' + 'A' * 9 + 'GTGA' + 'C' * 9 + 'GACA' + 'G' * 9 + 'T' * 8


def test_detect_bead_version_v1(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [V1_READ] * 100)
    assert wf.detect_bead_version(p) == 'v1'


def test_detect_bead_version_v1_gz(tmp_path):
    p = tmp_path / 'r1.fastq.gz'
    _write_fastq_gz(p, [V1_READ] * 100)
    assert wf.detect_bead_version(p) == 'v1'


def test_detect_bead_version_enhanced(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [ENH_READ] * 100)
    assert wf.detect_bead_version(p) == 'enhanced'


def test_detect_bead_version_enhanced_vb_stagger(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [ENH_READ_VB1] * 100)
    assert wf.detect_bead_version(p) == 'enhanced'


def test_detect_bead_version_unknown(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, ['ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'] * 100)
    assert wf.detect_bead_version(p) == 'unknown'


def test_detect_bead_version_empty(tmp_path):
    p = tmp_path / 'r1.fastq'
    p.write_text('')
    assert wf.detect_bead_version(p) == 'unknown'


def test_detect_bead_version_n_reads_limit(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [V1_READ] * 5 + ['ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'] * 1000)
    assert wf.detect_bead_version(p, n_reads=5) == 'v1'


## scan_r1_linkers histogram output (linker QC) -------------------------------

def test_scan_r1_linkers_v1_exact(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [V1_READ] * 100)
    stats = wf.scan_r1_linkers(p, n_reads=100)
    assert stats['class'] == 'v1'
    assert stats['n_reads'] == 100
    assert stats['v1_frac'] == 1.0
    assert stats['v1_errors'] == {0: 100}


def test_scan_r1_linkers_enhanced_exact(tmp_path):
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [ENH_READ] * 50)
    stats = wf.scan_r1_linkers(p, n_reads=50)
    assert stats['class'] == 'enhanced'
    assert stats['enh_frac'] == 1.0
    assert stats['enh_errors'][0] == 50


def test_scan_r1_linkers_v1_with_errors(tmp_path):
    ## introduce a single substitution in the first v1 linker (position 9)
    mutated = 'A' * 9 + 'TCTGGCCTGCGA' + 'C' * 9 + 'GGTAGCGGTGACA' + 'G' * 9 + 'T' * 8
    p = tmp_path / 'r1.fastq'
    _write_fastq(p, [V1_READ] * 80 + [mutated] * 20)
    stats = wf.scan_r1_linkers(p, n_reads=100)
    assert stats['n_reads'] == 100
    assert stats['v1_errors'].get(0) == 80
    assert stats['v1_errors'].get(1) == 20


def test_scan_r1_linkers_empty(tmp_path):
    p = tmp_path / 'r1.fastq'
    p.write_text('')
    stats = wf.scan_r1_linkers(p)
    assert stats['class'] == 'unknown'
    assert stats['n_reads'] == 0
    assert stats['v1_errors'] == {}
    assert stats['enh_errors'] == {}
