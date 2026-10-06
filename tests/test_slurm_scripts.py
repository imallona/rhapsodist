"""Checks on the job scripts under slurm/."""

import glob
import os
import re
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLURM_DIR = os.path.join(REPO_ROOT, 'slurm')
JOB_SCRIPTS = sorted(glob.glob(os.path.join(SLURM_DIR, '[0-9][0-9]_*.sh')))
REQUIRED_OPTIONS = ['job-name', 'time', 'cpus-per-task', 'mem-per-cpu', 'output']
## runs with benchmark files compared with each other
TIMED_STAGES = ['02', '03', '04', '05', '06', '07']

DIRECTIVE = re.compile(r'^#SBATCH --([a-z-]+)=(\S+)', re.MULTILINE)
CONFIG_PATH = re.compile(r'configs/[\w.-]+\.yaml')


def read(path):
    with open(path) as fh:
        return fh.read()


def sbatch_options(path):
    return dict(DIRECTIVE.findall(read(path)))


def stage(path):
    return os.path.basename(path)[:2]


def test_stages_are_numbered_from_zero():
    assert [stage(path) for path in JOB_SCRIPTS] == [
        f'{number:02d}' for number in range(len(JOB_SCRIPTS))]


@pytest.mark.parametrize('path', JOB_SCRIPTS + [os.path.join(SLURM_DIR, 'common.sh')],
                         ids=os.path.basename)
def test_script_parses(path):
    subprocess.run(['bash', '-n', path], check=True)


@pytest.mark.parametrize('path', JOB_SCRIPTS, ids=os.path.basename)
def test_script_requests_its_resources(path):
    options = sbatch_options(path)
    assert [name for name in REQUIRED_OPTIONS if name not in options] == []


def test_timed_runs_share_cpu_model_and_cores():
    timed = [sbatch_options(path) for path in JOB_SCRIPTS if stage(path) in TIMED_STAGES]
    assert len(timed) == len(TIMED_STAGES)
    assert len({options.get('constraint') for options in timed}) == 1
    assert timed[0].get('constraint')
    assert len({options['cpus-per-task'] for options in timed}) == 1


@pytest.mark.parametrize('path', JOB_SCRIPTS, ids=os.path.basename)
def test_named_configs_exist(path):
    for config in CONFIG_PATH.findall(read(path)):
        assert os.path.exists(os.path.join(REPO_ROOT, config)), config


def test_readme_lists_every_script():
    readme = read(os.path.join(REPO_ROOT, 'README.md'))
    missing = [os.path.basename(path) for path in JOB_SCRIPTS
               if os.path.basename(path) not in readme]
    assert missing == []
