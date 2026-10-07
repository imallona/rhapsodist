"""Dry runs of paper/Snakefile for the configs of the figure job, on a tree of
empty input files laid out as the workflow writes them."""

import os
import re
import shutil
import subprocess

import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURE_JOB = os.path.join(REPO_ROOT, 'slurm', '08_figures.sh')
DEFAULT_ALIGNERS = ['starsolo', 'kallisto', 'alevin']

pytestmark = pytest.mark.skipif(shutil.which('snakemake') is None, reason='needs snakemake')


def figure_configs():
    with open(FIGURE_JOB) as fh:
        return sorted(set(re.findall(r'configs/[\w.-]+\.yaml', fh.read())))


def touch(root, *parts):
    path = os.path.join(root, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, 'a').close()


def fake_inputs(root, config):
    wd = config['working_dir']
    samples = [s['name'] for s in config['samples']]
    aligners = DEFAULT_ALIGNERS + (['sbg'] if config.get('has_sbg') else [])
    runs = [wd] + list((config.get('linker_tolerance_runs') or {}).values())
    for sample in samples:
        for run in runs:
            for aligner in aligners:
                touch(root, run, aligner, sample, f'{sample}_{aligner}_sce.rds')
        touch(root, wd, 'linker_qc', f'{sample}_linker_errors.tsv')
        touch(root, wd, f'{sample}_biology_qc.rds')
        touch(root, wd, 'sampletags', sample, 'sampletag_demux.tsv.gz')
    touch(root, wd, 'run_info.tsv')


@pytest.mark.parametrize('config_path', figure_configs())
def test_paper_dry_run(config_path, tmp_path):
    with open(os.path.join(REPO_ROOT, config_path)) as fh:
        config = yaml.safe_load(fh)
    fake_inputs(tmp_path, config)
    result = subprocess.run(
        ['snakemake', '--snakefile', os.path.join(REPO_ROOT, 'paper', 'Snakefile'),
         '--configfile', os.path.join(REPO_ROOT, config_path),
         '--dry-run', '--nolock', '--cores', '1'],
        cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert re.search(r'^all\s+1$', result.stdout, re.M), result.stdout
