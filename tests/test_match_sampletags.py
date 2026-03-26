import gzip
import os
import tempfile

import pytest

from src.match_sampletags import best_match, load_tag_variable_regions

HUMAN_FA = os.path.join("workflow", "data", "sampletags", "human_sampletags.fa")
MOUSE_FA = os.path.join("workflow", "data", "sampletags", "mouse_sampletags.fa")

COMMON_PREFIX = "GTTGTCAAGATGCTACCGTTCAGAG"


def test_load_human_tags():
    tags = load_tag_variable_regions(HUMAN_FA)
    assert len(tags) == 12
    assert all(isinstance(v, str) and len(v) > 0 for v in tags.values())
    assert "human_sampletag_1" in tags


def test_load_mouse_tags():
    tags = load_tag_variable_regions(MOUSE_FA)
    assert len(tags) == 12
    assert "mouse_sampletag_1" in tags


def test_variable_regions_are_distinct():
    tags = load_tag_variable_regions(HUMAN_FA)
    regions = list(tags.values())
    assert len(set(regions)) == len(regions)


def test_exact_match():
    tags = load_tag_variable_regions(HUMAN_FA)
    first_name = next(iter(tags))
    first_var = tags[first_name]
    tag, n_mm = best_match(first_var, tags, max_mismatch=2)
    assert tag == first_name
    assert n_mm == 0


def test_one_mismatch():
    tags = load_tag_variable_regions(HUMAN_FA)
    first_name = next(iter(tags))
    first_var = list(tags[first_name])
    first_var[0] = "N"
    mutated = "".join(first_var)
    tag, n_mm = best_match(mutated, tags, max_mismatch=2)
    assert tag == first_name
    assert n_mm == 1


def test_too_many_mismatches_returns_none():
    tags = load_tag_variable_regions(HUMAN_FA)
    tag, n_mm = best_match("N" * 45, tags, max_mismatch=2)
    assert tag is None
    assert n_mm is None


def test_main_writes_output(tmp_path):
    tags = load_tag_variable_regions(HUMAN_FA)
    first_name = next(iter(tags))
    seq = COMMON_PREFIX + tags[first_name]

    fq = tmp_path / "test.fq.gz"
    out = tmp_path / "out.tsv.gz"

    with gzip.open(fq, "wt") as fh:
        fh.write(f"@ACGT__TTTT\n{seq}\n+\n{'I' * len(seq)}\n")
        fh.write(f"@GGGG__CCCC\nACGTACGT\n+\nIIIIIIII\n")

    from src.match_sampletags import main
    import sys

    argv_backup = sys.argv
    sys.argv = ["match_sampletags", "--fa", HUMAN_FA,
                "--fastq", str(fq), "--out", str(out), "--max-mismatch", "2"]
    try:
        main()
    except SystemExit:
        pass
    finally:
        sys.argv = argv_backup

    with gzip.open(out, "rt") as fh:
        lines = fh.readlines()

    assert len(lines) == 1
    fields = lines[0].strip().split("\t")
    assert fields[0] == "ACGT"
    assert fields[1] == "TTTT"
    assert fields[2] == first_name
    assert fields[3] == "0"
