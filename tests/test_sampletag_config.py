import os.path
import types

import pytest

from workflow.src import workflow_functions as wf

# workflow_functions.py is normally include()d into the Snakefile namespace, so the
# globals it reads (op, config, workflow) are injected here for unit testing.
wf.op = os.path
wf.workflow = types.SimpleNamespace(basedir="workflow")


def _set_config(sampletags=None, species="human"):
    uses = {"species": species, "use_sampletags": "yes"}
    if sampletags is not None:
        uses["sampletags"] = sampletags
    wf.config = {"skip_sampletags": False, "samples": [{"name": "s1", "uses": uses}]}


def test_list_form_labels_default_to_tag_name():
    _set_config([1, 2])
    tags = wf.get_sampletags_by_name("s1")
    assert tags == {"human_sampletag_1": "human_sampletag_1",
                    "human_sampletag_2": "human_sampletag_2"}


def test_mapping_form_renames_labels():
    _set_config({1: "donorA", 2: "donorB"})
    tags = wf.get_sampletags_by_name("s1")
    assert tags == {"human_sampletag_1": "donorA", "human_sampletag_2": "donorB"}


def test_full_tag_name_is_accepted():
    _set_config(["human_sampletag_3"])
    assert wf.get_sampletag_labels_by_name("s1") == ["human_sampletag_3"]


def test_absent_field_returns_none():
    _set_config(None)
    assert wf.get_sampletags_by_name("s1") is None
    assert wf.get_sampletag_labels_by_name("s1") is None


def test_unknown_tag_raises():
    _set_config([999])
    with pytest.raises(ValueError):
        wf.get_sampletags_by_name("s1")


def test_sampletags_without_species_raises():
    _set_config([1], species=None)
    # species None is rejected by get_species_by_name before tag resolution
    with pytest.raises(ValueError):
        wf.get_sampletags_by_name("s1")
