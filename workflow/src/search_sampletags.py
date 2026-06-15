#!/usr/bin/env python
##
## Alignment-free sampletag calling by direct sequence search.
##
## Used only when starsolo is not among the configured aligners. When starsolo
## runs, sampletags are called by the published starsolo mode instead (see the
## extract_unmapped_startsolo_wta_tagged_fastqs ... count_sampletags rules), so
## that mode and its figures stay reproducible.
##
## A BD sampletag read carries a fixed prefix shared by all tags followed by a
## tag-specific region. We anchor on the prefix (allowing a few mismatches),
## assign the tag region to the closest reference tag, and recover the cell
## barcode and UMI from the paired standardized cb/umi read. Cell barcode
## segments are corrected against the BD whitelists with the same one-mismatch
## tolerance starsolo applies, so the barcodes match the cells the aligner calls.
##
## Output is the same four-column table the starsolo mode emits via
## match_sampletags (cb, umi, sampletag, n_mismatches), so the downstream demux
## and report are shared and unchanged.
##
## GPLv3

import argparse
import gzip

try:
    ## when run as a Snakemake script, the sibling module is on sys.path
    from match_sampletags import common_prefix, load_tag_variable_regions, best_match
except ImportError:
    ## when imported as a package (pytest)
    from workflow.src.match_sampletags import common_prefix, load_tag_variable_regions, best_match


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def build_whitelist_corrector(barcodes):
    """Map each whitelist barcode and each of its one-mismatch neighbors to the
    canonical barcode. Neighbors reachable from two different barcodes are
    ambiguous and dropped, so a lookup only succeeds when correction is unique."""
    corrector = {}
    ambiguous = set()
    bases = 'ACGT'
    for bc in barcodes:
        corrector[bc] = bc
        for i, ref in enumerate(bc):
            for b in bases:
                if b == ref:
                    continue
                neighbor = bc[:i] + b + bc[i + 1:]
                if neighbor in corrector and corrector[neighbor] != bc:
                    ambiguous.add(neighbor)
                else:
                    corrector[neighbor] = bc
    for n in ambiguous:
        del corrector[n]
    return corrector


def correct_barcode(cb, correctors, segment_len):
    """Correct a concatenated cell barcode segment by segment against the
    per-segment whitelist correctors. Return the corrected barcode, or None when
    any segment cannot be resolved to a unique whitelist entry."""
    if len(cb) != segment_len * len(correctors):
        return None
    out = []
    for i, corrector in enumerate(correctors):
        seg = cb[i * segment_len:(i + 1) * segment_len]
        canonical = corrector.get(seg)
        if canonical is None:
            return None
        out.append(canonical)
    return ''.join(out)


def _open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def _fastq_seqs(path):
    """Yield the sequence line of each record in a fastq."""
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip().upper()


def search(r1_path, r2_path, tags_path, whitelist_paths, out_path,
           segment_len=9, umi_len=8, max_anchor_dist=2, max_tag_dist=2):
    """Stream the paired standardized reads, call sampletags by search, and write
    the cb/umi/sampletag/n_mismatches table. r1 carries the cell barcode (segments
    of segment_len) followed by the UMI; r2 carries the tag read."""
    tag_vars = load_tag_variable_regions(tags_path)
    correctors = [build_whitelist_corrector(_read_lines(p)) for p in whitelist_paths]
    cb_len = segment_len * len(correctors)
    anchor_len = len(common_prefix)

    n_total = n_anchored = n_assigned = n_written = 0
    with gzip.open(out_path, 'wt') as out:
        for cbumi, r2 in zip(_fastq_seqs(r1_path), _fastq_seqs(r2_path)):
            n_total += 1
            if len(r2) < anchor_len:
                continue
            if hamming(r2[:anchor_len], common_prefix) > max_anchor_dist:
                continue
            n_anchored += 1
            tag, n_mm = best_match(r2[anchor_len:], tag_vars, max_tag_dist)
            if tag is None:
                continue
            n_assigned += 1
            cb = correct_barcode(cbumi[:cb_len], correctors, segment_len)
            if cb is None:
                continue
            umi = cbumi[cb_len:cb_len + umi_len]
            out.write(f"{cb}\t{umi}\t{tag}\t{n_mm}\n")
            n_written += 1

    return {'reads': n_total, 'anchored': n_anchored,
            'assigned': n_assigned, 'written': n_written}


def _read_lines(path):
    with open(path) as fh:
        return [line.strip().upper() for line in fh if line.strip()]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--r1', required=True, help='standardized cb/umi fastq')
    p.add_argument('--r2', required=True, help='standardized cdna fastq')
    p.add_argument('--tags', required=True, help='species sampletag fasta')
    p.add_argument('--whitelists', required=True, nargs='+',
                   help='per-segment BD_CLS whitelist files, in barcode order')
    p.add_argument('--out', required=True, help='output counts tsv.gz')
    p.add_argument('--segment-len', type=int, default=9)
    p.add_argument('--umi-len', type=int, default=8)
    p.add_argument('--max-anchor-dist', type=int, default=2)
    p.add_argument('--max-tag-dist', type=int, default=2)
    p.add_argument('--log', default=None)
    args = p.parse_args()

    stats = search(args.r1, args.r2, args.tags, args.whitelists, args.out,
                   segment_len=args.segment_len, umi_len=args.umi_len,
                   max_anchor_dist=args.max_anchor_dist, max_tag_dist=args.max_tag_dist)

    summary = (
        f"reads scanned:               {stats['reads']:>10}\n"
        f"anchor matched:              {stats['anchored']:>10}\n"
        f"tag assigned:                {stats['assigned']:>10}\n"
        f"written (barcode corrected): {stats['written']:>10}\n"
    )
    if args.log:
        with open(args.log, 'w') as fh:
            fh.write(summary)
    print(summary, end='')


if __name__ == '__main__':
    main()
