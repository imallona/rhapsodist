import gzip
import os

import pytest

from workflow.src.search_sampletags import (
    build_whitelist_corrector,
    correct_barcode,
    hamming,
    search,
)
from workflow.src.match_sampletags import common_prefix, load_tag_variable_regions

HUMAN_FA = os.path.join("workflow", "data", "sampletags", "human_sampletags.fa")


def test_hamming():
    assert hamming("AAAA", "AAAA") == 0
    assert hamming("AAAA", "AAAT") == 1
    assert hamming("AAAA", "TTTT") == 4


def test_whitelist_corrector_exact_and_one_mismatch():
    corrector = build_whitelist_corrector(["AAAAAAAAA", "CCCCCCCCC"])
    assert corrector["AAAAAAAAA"] == "AAAAAAAAA"
    # one substitution away from the first barcode
    assert corrector["TAAAAAAAA"] == "AAAAAAAAA"
    # two substitutions away resolves to nothing
    assert "TTAAAAAAA" not in corrector


def test_whitelist_corrector_drops_ambiguous_neighbors():
    # AAAAAAAAA and AAAAAAAAT are one apart; the neighbor that sits between them
    # is reachable from both and must be dropped
    corrector = build_whitelist_corrector(["AAAAAAAAA", "TAAAAAAAT"])
    # AAAAAAAAT is 1 from AAAAAAAAA (last base) and 1 from TAAAAAAAT (first base)
    assert "AAAAAAAAT" not in corrector


def test_correct_barcode_three_segments():
    correctors = [
        build_whitelist_corrector(["AAAAAAAAA"]),
        build_whitelist_corrector(["CCCCCCCCC"]),
        build_whitelist_corrector(["GGGGGGGGG"]),
    ]
    cb = "AAAAAAAAA" + "CCCCCCCCC" + "GGGGGGGGG"
    assert correct_barcode(cb, correctors, 9) == cb
    # one error in the middle segment is corrected
    cb_err = "AAAAAAAAA" + "CCCCCTCCC" + "GGGGGGGGG"
    assert correct_barcode(cb_err, correctors, 9) == cb


def test_correct_barcode_rejects_unresolvable():
    correctors = [build_whitelist_corrector(["AAAAAAAAA"])]
    assert correct_barcode("TTTTTTTTT", correctors, 9) is None
    # wrong total length
    assert correct_barcode("AAAA", correctors, 9) is None


def _write_fastq(path, records):
    with gzip.open(path, "wt") as fh:
        for name, seq in records:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def test_search_end_to_end(tmp_path):
    tag_vars = load_tag_variable_regions(HUMAN_FA)
    tag1 = "human_sampletag_1"
    tag1_region = tag_vars[tag1]

    cb_seg = "AAAAAAAAA"
    cb = cb_seg * 3
    umi = "ACGTACGT"
    r1 = tmp_path / "r1.fq.gz"
    r2 = tmp_path / "r2.fq.gz"

    # read 1: a clean sampletag read for tag1
    # read 2: a non-tag read (no anchor) that must be dropped
    # read 3: tag1 read whose barcode has one correctable error
    _write_fastq(r1, [
        ("r1", cb + umi),
        ("r2", cb + umi),
        ("r3", "TAAAAAAAA" + cb_seg + cb_seg + umi),
    ])
    _write_fastq(r2, [
        ("r1", common_prefix + tag1_region),
        ("r2", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
        ("r3", common_prefix + tag1_region),
    ])

    whitelist = tmp_path / "BD_CLS1.txt"
    whitelist.write_text(cb_seg + "\n")
    out = tmp_path / "counts.tsv.gz"

    stats = search(str(r1), str(r2), HUMAN_FA,
                   [str(whitelist), str(whitelist), str(whitelist)], str(out))

    assert stats["reads"] == 3
    assert stats["anchored"] == 2
    assert stats["assigned"] == 2
    assert stats["written"] == 2

    with gzip.open(out, "rt") as fh:
        rows = [line.strip().split("\t") for line in fh]
    assert len(rows) == 2
    for r in rows:
        assert r[0] == cb           # both barcodes corrected to the canonical
        assert r[1] == umi
        assert r[2] == tag1
