#!/usr/bin/env python3
"""
Identify the knee-point on the alevin barcode rank plot (DeduplicatedReads vs rank)
and write the barcodes above the knee to a file.

For simulated data, pass --sim_barcodes to copy the known cell list directly.
"""

import argparse
import math
import shutil
import sys


def find_knee_threshold(counts_sorted_desc):
    """
    Kneedle algorithm on the log10-log10 barcode rank plot.
    Returns the count value at the knee (inclusive lower bound).
    counts_sorted_desc: list of counts in descending order, already filtered to > 0.
    """
    counts = [c for c in counts_sorted_desc if c > 0]
    if len(counts) < 3:
        return counts[-1] if counts else 0.0

    log_r = [math.log10(i + 1) for i in range(len(counts))]
    log_c = [math.log10(c) for c in counts]

    r_min, r_max = log_r[0], log_r[-1]
    c_min, c_max = log_c[-1], log_c[0]

    if r_max == r_min or c_max == c_min:
        return float(counts[-1])

    x = [(r - r_min) / (r_max - r_min) for r in log_r]
    y = [(c - c_min) / (c_max - c_min) for c in log_c]

    # perpendicular distance from the line (0,1) -> (1,0): x + y = 1
    dists = [abs(xi + yi - 1) / math.sqrt(2) for xi, yi in zip(x, y)]
    knee_idx = dists.index(max(dists))
    return float(counts[knee_idx])


def main():
    parser = argparse.ArgumentParser(
        description='Knee-point barcode filter from alevin featureDump.txt')
    parser.add_argument('--feature_dump', required=True,
                        help='Path to alevin featureDump.txt')
    parser.add_argument('--output', required=True,
                        help='Output path for filtered barcodes (one per line)')
    parser.add_argument('--sim_barcodes', default='',
                        help='Simulation mode: copy this file to output, skipping knee detection')
    args = parser.parse_args()

    if args.sim_barcodes:
        shutil.copy(args.sim_barcodes, args.output)
        print(f'simulation mode: copied {args.sim_barcodes} to {args.output}',
              file=sys.stderr)
        return

    barcodes = []
    dedup_reads = []
    with open(args.feature_dump) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            barcodes.append(parts[0])
            dedup_reads.append(float(parts[3]))

    counts_sorted = sorted(dedup_reads, reverse=True)
    threshold = find_knee_threshold(counts_sorted)

    kept = [b for b, c in zip(barcodes, dedup_reads) if c >= threshold]
    with open(args.output, 'w') as out:
        for b in kept:
            out.write(b + '\n')

    print(
        f'knee threshold: {threshold:.1f} deduplicated reads, '
        f'kept {len(kept)} / {len(barcodes)} barcodes',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
