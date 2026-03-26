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
        'whitelist': '96x3',
        'species': 'human',
    },
}

SAMPLE_B = {
    'name': 'sampleB',
    'uses': {
        'cb_umi_fq': '/data/sampleB_R1.fq.gz',
        'cdna_fq': '/data/sampleB_R2.fq.gz',
        'whitelist': '384x3',
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


def test_get_barcode_whitelist_96x3():
    assert wf.get_barcode_whitelist_by_name('sampleA') == '96x3'


def test_get_barcode_whitelist_384x3():
    assert wf.get_barcode_whitelist_by_name('sampleB') == '384x3'


def test_get_species_human():
    assert wf.get_species_by_name('sampleA') == 'human'


def test_get_species_mouse():
    assert wf.get_species_by_name('sampleB') == 'mouse'


def test_get_species_invalid_raises():
    wf.config['samples'][0]['uses']['species'] = 'zebrafish'
    with pytest.raises(ValueError):
        wf.get_species_by_name('sampleA')


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
