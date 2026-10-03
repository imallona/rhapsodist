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


def _set_library_config(**uses):
    base = {"species": "mouse", "use_sampletags": "yes", "allowedlist": 384,
            "diversity_insets": "yes", "cb_umi_fq": "/d/wta_R1.fq.gz",
            "cdna_fq": "/d/wta_R2.fq.gz", "downsample": 10}
    base.update(uses)
    wf.config = {"skip_sampletags": False, "working_dir": "/wd", "aligner": ["starsolo", "alevin"],
                 "samples": [{"name": "s1", "uses": base}]}


def test_sample_without_tag_library_reads_tags_from_itself():
    _set_library_config()
    assert not wf.has_sampletag_library("s1")
    assert wf.sampletag_reads_name("s1") == "s1"
    assert wf.sampletag_counts_by_name("s1").endswith("sampletag_counts.tsv.gz")


def test_tag_library_fastqs_become_a_library_entry():
    _set_library_config(sampletag_cb_umi_fq="/d/st_R1.fq.gz", sampletag_cdna_fq="/d/st_R2.fq.gz")
    assert wf.has_sampletag_library("s1")
    assert wf.sampletag_reads_name("s1") == "s1_sampletags"
    assert wf.get_cbumi_by_name("s1_sampletags") == "/d/st_R1.fq.gz"
    assert wf.get_cdna_by_name("s1_sampletags") == "/d/st_R2.fq.gz"
    assert wf.get_declared_allowedlist("s1_sampletags") == wf.get_declared_allowedlist("s1")


def test_tag_library_is_not_downsampled():
    _set_library_config(sampletag_cb_umi_fq="/d/st_R1.fq.gz", sampletag_cdna_fq="/d/st_R2.fq.gz")
    assert wf.get_downsample_fraction("s1") == 0.1
    assert wf.get_downsample_fraction("s1_sampletags") == 1.0


def test_tag_library_is_not_a_sample():
    _set_library_config(sampletag_cb_umi_fq="/d/st_R1.fq.gz", sampletag_cdna_fq="/d/st_R2.fq.gz")
    assert wf.get_sample_names() == ["s1"]
    assert wf.samples_with_sampletags() == ["s1"]


def test_tag_library_uses_search_counts_even_with_starsolo():
    _set_library_config(sampletag_cb_umi_fq="/d/st_R1.fq.gz", sampletag_cdna_fq="/d/st_R2.fq.gz")
    assert wf.get_sampletag_method() == "starsolo"
    assert wf.sampletag_counts_by_name("s1").endswith("sampletag_counts_search.tsv.gz")


def test_tag_library_from_sra():
    _set_library_config(sampletag_sra_run="SRR0000001")
    assert wf.has_sampletag_library("s1")
    assert wf.get_cbumi_by_name("s1_sampletags") == \
        "/wd/data/fastq/sra/s1_sampletags/s1_sampletags_R1.fastq.gz"


def test_tag_library_with_one_mate_raises():
    _set_library_config(sampletag_cb_umi_fq="/d/st_R1.fq.gz")
    with pytest.raises(ValueError, match="given together"):
        wf.validate_sampletag_library("s1")


def test_tag_library_with_fastqs_and_sra_raises():
    _set_library_config(sampletag_cb_umi_fq="/d/st_R1.fq.gz", sampletag_cdna_fq="/d/st_R2.fq.gz",
                        sampletag_sra_run="SRR0000001")
    with pytest.raises(ValueError, match="not both"):
        wf.validate_sampletag_library("s1")


def test_tag_library_name_clash_raises():
    _set_library_config(sampletag_sra_run="SRR0000001")
    wf.config["samples"].append({"name": "s1_sampletags", "uses": {"cb_umi_fq": "a", "cdna_fq": "b"}})
    with pytest.raises(ValueError, match="clashes"):
        wf.validate_sampletag_library("s1")
