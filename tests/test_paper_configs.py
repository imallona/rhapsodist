"""The comparison keys of the configs point at runs and files that exist."""

import glob
import os

import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIGS = sorted(glob.glob(os.path.join(REPO_ROOT, 'configs', '*_config.yaml')))
KNOWN_WORKFLOWS = {'universc', 'zumis', 'openpipelines'}


def load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def configs_with(key):
    return [path for path in CONFIGS if key in load(path)]


@pytest.mark.parametrize('path', configs_with('linker_tolerance_runs'), ids=os.path.basename)
def test_linker_runs_match_their_configs(path):
    stem = os.path.basename(path).replace('_config.yaml', '')
    for errors, working_dir in load(path)['linker_tolerance_runs'].items():
        run = load(os.path.join(REPO_ROOT, 'configs', f'{stem}_linker{errors}_config.yaml'))
        assert run['working_dir'] == working_dir
        assert run['cb_umi_max_errors'] == errors


@pytest.mark.parametrize('path', configs_with('other_workflows'), ids=os.path.basename)
def test_other_workflows_are_known(path):
    assert set(load(path)['other_workflows']) <= KNOWN_WORKFLOWS


@pytest.mark.parametrize('path', configs_with('deposited_sampletag_files'), ids=os.path.basename)
def test_deposited_files_cover_the_sample_tags(path):
    config = load(path)
    deposited = config['deposited_sampletag_files']
    for sample in config['samples']:
        uses = sample['uses']
        expected = {f"{uses['species']}_sampletag_{tag}" for tag in uses['sampletags']}
        assert set(deposited) == expected
    assert all(url.startswith('https://') for url in deposited.values())
