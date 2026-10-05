import subprocess
import sys
from types import SimpleNamespace

import pytest

from rhapsodist import cli


def snakemake_command(monkeypatch, *arguments):
    calls = []
    monkeypatch.setattr(sys, 'argv', ['rhapsodist', *arguments])
    monkeypatch.setattr(subprocess, 'run',
                        lambda cmd: calls.append(cmd) or SimpleNamespace(returncode=0))
    with pytest.raises(SystemExit) as exit_info:
        cli.main()
    assert exit_info.value.code == 0
    return calls[0]


def test_config_overrides_are_forwarded_to_snakemake(monkeypatch):
    cmd = snakemake_command(monkeypatch, '--configfile', 'a.yaml', '--config', 'aligner=x')
    assert cmd[cmd.index('--configfile') + 1] == 'a.yaml'
    assert cmd[-2:] == ['--config', 'aligner=x']


def test_dry_run_and_cores(monkeypatch):
    cmd = snakemake_command(monkeypatch, '--configfile', 'a.yaml', '--cores', '3', '-n')
    assert cmd[cmd.index('--cores') + 1] == '3'
    assert '--dry-run' in cmd
