import os.path
import types

import pytest

from workflow.src import workflow_functions as wf

# workflow_functions.py is normally include()d into the Snakefile namespace, so the
# globals it reads (op, config, workflow) are injected here for unit testing.
wf.op = os.path
wf.workflow = types.SimpleNamespace(basedir="workflow")


def _set_config(output_format=None, backend=None, working_dir="output",
                aligners=("starsolo", "alevin"), samples=("s1",)):
    cfg = {"working_dir": working_dir, "aligner": list(aligners),
           "samples": [{"name": s, "uses": {}} for s in samples]}
    if output_format is not None:
        cfg["output_format"] = output_format
    if backend is not None:
        cfg["sampletag_split_backend"] = backend
    wf.config = cfg


def test_output_format_defaults_to_sce():
    _set_config()
    assert wf.get_output_format() == "sce"
    assert wf.wants_h5ad() is False


@pytest.mark.parametrize("fmt,wants", [("sce", False), ("h5ad", True), ("both", True)])
def test_output_format_values(fmt, wants):
    _set_config(output_format=fmt)
    assert wf.get_output_format() == fmt
    assert wf.wants_h5ad() is wants


def test_output_format_is_case_insensitive():
    _set_config(output_format="BOTH")
    assert wf.get_output_format() == "both"


def test_invalid_output_format_raises():
    _set_config(output_format="loom")
    with pytest.raises(ValueError):
        wf.get_output_format()


def test_split_backend_defaults_to_memory():
    _set_config()
    assert wf.get_sampletag_split_backend() == "memory"


@pytest.mark.parametrize("backend", ["memory", "delayed"])
def test_split_backend_values(backend):
    _set_config(backend=backend)
    assert wf.get_sampletag_split_backend() == backend


def test_invalid_split_backend_raises():
    _set_config(backend="streaming")
    with pytest.raises(ValueError):
        wf.get_sampletag_split_backend()


def test_split_flags_compose_backend_and_format():
    _set_config(output_format="both", backend="delayed")
    assert wf.sampletag_split_flags() == "--backend delayed --output_format both"


def test_h5ad_targets_empty_when_sce_only():
    _set_config(output_format="sce")
    assert wf.h5ad_sce_targets() == []


def test_h5ad_targets_cover_every_aligner_and_sample():
    _set_config(output_format="both", working_dir="out",
                aligners=("starsolo", "alevin"), samples=("a", "b"))
    targets = set(wf.h5ad_sce_targets())
    assert targets == {
        "out/starsolo/a/a_starsolo.h5ad",
        "out/alevin/a/a_alevin.h5ad",
        "out/starsolo/b/b_starsolo.h5ad",
        "out/alevin/b/b_alevin.h5ad",
    }
