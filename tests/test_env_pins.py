import glob
import os
import re

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENV_DIR = os.path.join(REPO_ROOT, 'workflow', 'envs')
PINNED_ENVS = sorted(
    path for path in glob.glob(os.path.join(ENV_DIR, '*.yaml'))
    if os.path.exists(path[:-len('.yaml')] + '.linux-64.pin.txt')
)

DEPENDENCY = re.compile(r'^\s*-\s*(?:[\w-]+::)?([A-Za-z0-9_.-]+)')


def yaml_packages(path):
    """Conda package names under dependencies, without channel and version."""
    packages = []
    in_dependencies = False
    with open(path) as fh:
        for line in fh:
            if line.startswith('dependencies:'):
                in_dependencies = True
                continue
            match = DEPENDENCY.match(line.split('#')[0])
            if in_dependencies and match:
                packages.append(match.group(1).lower())
    return packages


def pinned_packages(path):
    """Package names of an explicit conda pin file, from its package URLs."""
    names = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('http'):
                filename = line.strip().rsplit('/', 1)[1]
                names.add(filename.rsplit('-', 2)[0].lower())
    return names


def test_pinned_environments_are_found():
    assert len(PINNED_ENVS) >= 8


@pytest.mark.parametrize('yaml_path', PINNED_ENVS, ids=os.path.basename)
def test_every_yaml_package_is_pinned(yaml_path):
    pinned = pinned_packages(yaml_path[:-len('.yaml')] + '.linux-64.pin.txt')
    missing = [name for name in yaml_packages(yaml_path) if name not in pinned]
    assert missing == []
