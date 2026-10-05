#!/usr/bin/env python3
"""Cell barcode matches against the whitelist, as each aligner reports them.

The aligners count different things: STARsolo counts reads, bustools counts BUS
records, alevin-fry counts distinct barcodes it corrected, and salmon alevin
without --sketch counts only the reads it discarded. A count a tool does not
report is left empty.
"""

import argparse
import json
import re

COLUMNS = ['aligner', 'unit', 'exact', 'corrected', 'unmatched']
STAR_CORRECTED = ('yesOneWLmatchWithMM', 'yesMultWLmatchWithMM')
STAR_UNMATCHED = ('noNoWLmatch', 'noTooManyMM', 'noTooManyWLmatches')
BUSTOOLS_COUNT = re.compile(r'^(In on-list|Corrected|Uncorrected)\s*=\s*(\d+)\s*$')
FRY_CORRECTED = re.compile(r'total number of distinct corrected barcodes\s*:\s*([\d,]+)')


def starsolo_counts(lines):
    """From Solo.out/Barcodes.stats: one name and one read count per line."""
    counts = {}
    for line in lines:
        fields = line.split()
        if len(fields) == 2:
            counts[fields[0]] = int(fields[1])
    return {'aligner': 'starsolo', 'unit': 'reads',
            'exact': counts['yesWLmatchExact'],
            'corrected': sum(counts[name] for name in STAR_CORRECTED),
            'unmatched': sum(counts[name] for name in STAR_UNMATCHED)}


def bustools_counts(lines):
    """From the log of bustools correct."""
    counts = {}
    for line in lines:
        match = BUSTOOLS_COUNT.match(line.strip())
        if match:
            counts[match[1]] = int(match[2])
    return {'aligner': 'kallisto', 'unit': 'BUS records',
            'exact': counts['In on-list'],
            'corrected': counts['Corrected'],
            'unmatched': counts['Uncorrected']}


def alevin_fry_counts(lines):
    """From the log of alevin-fry generate-permit-list."""
    for line in lines:
        match = FRY_CORRECTED.search(line)
        if match:
            return {'aligner': 'alevin', 'unit': 'distinct barcodes', 'exact': '',
                    'corrected': int(match[1].replace(',', '')), 'unmatched': ''}
    raise ValueError('no corrected barcode count in the alevin-fry log')


def salmon_alevin_counts(meta):
    """From aux_info/alevin_meta_info.json of salmon alevin without --sketch."""
    return {'aligner': 'alevin', 'unit': 'reads', 'exact': '', 'corrected': '',
            'unmatched': meta['noisy_cb_reads']}


def write_counts(rows, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(COLUMNS) + '\n')
        for row in rows:
            fh.write('\t'.join(str(row[c]) for c in COLUMNS) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--starsolo_stats', help='Solo.out/Barcodes.stats')
    ap.add_argument('--bustools_log', help='log of bustools correct')
    ap.add_argument('--alevin_fry_log', help='log of alevin-fry generate-permit-list')
    ap.add_argument('--alevin_meta', help='aux_info/alevin_meta_info.json')
    ap.add_argument('--out', required=True, help='output TSV')
    args = ap.parse_args()

    rows = []
    for path, parse in ((args.starsolo_stats, starsolo_counts),
                        (args.bustools_log, bustools_counts),
                        (args.alevin_fry_log, alevin_fry_counts)):
        if path:
            with open(path, errors='replace') as fh:
                rows.append(parse(fh))
    if args.alevin_meta:
        with open(args.alevin_meta) as fh:
            rows.append(salmon_alevin_counts(json.load(fh)))
    write_counts(rows, args.out)


if __name__ == '__main__':
    main()
