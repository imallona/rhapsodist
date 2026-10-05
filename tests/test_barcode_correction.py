import json

import pytest

from workflow.src.barcode_correction import (alevin_fry_counts, bustools_counts, main,
                                             salmon_alevin_counts, starsolo_counts)

## Solo.out/Barcodes.stats of STAR 2.7.11b, 10 percent of SRR24978231.
STAR_STATS = """                                       noNoAdapter              0
                                           noNoUMI              0
                                            noNoCB              0
                                           noNinCB           3291
                                          noNinUMI            554
                                  noUMIhomopolymer         134074
                                       noNoWLmatch        5391116
                                       noTooManyMM              0
                                noTooManyWLmatches              0
                                   yesWLmatchExact       42049264
                               yesOneWLmatchWithMM         740748
                              yesMultWLmatchWithMM        2275402
"""

BUSTOOLS_LOG = """Found 391958 barcodes in the on-list
Processed 38543851 BUS records
In on-list = 32879850
Corrected    = 3710241
Uncorrected  = 1953760
"""

FRY_LOG = """2026-10-03 04:51:11 INFO observed 41,893,907 reads (41,154,653 orientation consistent) in 8,381 chunks
2026-10-03 04:51:27 INFO total number of distinct corrected barcodes : 2,672,696
2026-10-03 04:51:27 INFO filter_type = Filtered
"""


def test_starsolo_adds_both_kinds_of_corrected_reads():
    row = starsolo_counts(STAR_STATS.splitlines())
    assert row == {'aligner': 'starsolo', 'unit': 'reads', 'exact': 42049264,
                   'corrected': 740748 + 2275402, 'unmatched': 5391116}


def test_bustools_counts_bus_records():
    row = bustools_counts(BUSTOOLS_LOG.splitlines())
    assert row == {'aligner': 'kallisto', 'unit': 'BUS records', 'exact': 32879850,
                   'corrected': 3710241, 'unmatched': 1953760}
    assert row['exact'] + row['corrected'] + row['unmatched'] == 38543851


def test_alevin_fry_reports_distinct_barcodes_only():
    row = alevin_fry_counts(FRY_LOG.splitlines())
    assert row == {'aligner': 'alevin', 'unit': 'distinct barcodes', 'exact': '',
                   'corrected': 2672696, 'unmatched': ''}


def test_alevin_fry_log_without_the_count_is_an_error():
    with pytest.raises(ValueError):
        alevin_fry_counts(['INFO filter_type = Filtered'])


def test_salmon_alevin_reports_discarded_reads_only():
    row = salmon_alevin_counts({'total_reads': 404743, 'noisy_cb_reads': 274481})
    assert row == {'aligner': 'alevin', 'unit': 'reads', 'exact': '', 'corrected': '',
                   'unmatched': 274481}


def test_table_has_one_row_per_given_input(tmp_path, monkeypatch):
    stats = tmp_path / 'Barcodes.stats'
    stats.write_text(STAR_STATS)
    meta = tmp_path / 'alevin_meta_info.json'
    meta.write_text(json.dumps({'noisy_cb_reads': 7}))
    out = tmp_path / 'counts.tsv'
    monkeypatch.setattr('sys.argv', ['barcode_correction.py', '--starsolo_stats', str(stats),
                                     '--alevin_meta', str(meta), '--out', str(out)])
    main()
    assert out.read_text().splitlines() == [
        'aligner\tunit\texact\tcorrected\tunmatched',
        'starsolo\treads\t42049264\t3016150\t5391116',
        'alevin\treads\t\t\t7',
    ]
