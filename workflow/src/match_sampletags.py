#!/usr/bin/env python3
##
## fast sampletag demultiplexing by direct hamming-distance matching.
##
## reads a cb__umi-labelled fastq (from extract_unmapped_startsolo_wta_tagged_fastqs),
## filters reads starting with the shared sampletag prefix,
## assigns each read to the closest sampletag by hamming distance on the variable region,
## writes a tsv.gz: cb, umi, sampletag, n_mismatches.
##
## izaskun mallona

import argparse
import gzip
import sys

common_prefix = "GTTGTCAAGATGCTACCGTTCAGAG"


def load_tag_variable_regions(fa_path):
    """return {tag_name: variable_suffix} for each entry in the fasta."""
    tags = {}
    name = None
    with open(fa_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif name is not None:
                seq = line.upper()
                if not seq.startswith(common_prefix):
                    sys.exit(
                        f"tag {name!r} does not start with the expected common prefix.\n"
                        f"  expected: {common_prefix}\n"
                        f"  got:      {seq[:len(common_prefix)]}"
                    )
                tags[name] = seq[len(common_prefix):]
                name = None
    return tags


def best_match(variable_seq, tag_vars, max_mismatch):
    """return (tag_name, n_mismatches), or (None, None) if no tag is within threshold."""
    best_tag = None
    best_dist = max_mismatch + 1
    for tag, ref in tag_vars.items():
        n = min(len(variable_seq), len(ref))
        dist = sum(a != b for a, b in zip(variable_seq[:n], ref[:n]))
        if dist < best_dist:
            best_dist = dist
            best_tag = tag
    return (best_tag, best_dist) if best_tag is not None else (None, None)


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def main():
    ap = argparse.ArgumentParser(
        description="assign sampletag reads by hamming-distance matching on the variable region."
    )
    ap.add_argument("--fa", required=True, help="sampletag fasta (species-specific)")
    ap.add_argument("--fastq", required=True, help="fastq(.gz) with @cb__umi headers")
    ap.add_argument("--out", required=True, help="output tsv.gz: cb, umi, sampletag, n_mismatches")
    ap.add_argument("--max-mismatch", type=int, default=2, metavar="n",
                    help="maximum hamming distance allowed in the variable region (default: 2)")
    args = ap.parse_args()

    tag_vars = load_tag_variable_regions(args.fa)
    if not tag_vars:
        sys.exit("no tags loaded from fasta - check the file path.")

    n_read = n_prefix = n_matched = 0

    with open_maybe_gz(args.fastq) as fq, gzip.open(args.out, "wt") as out:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline().rstrip("\n")
            fq.readline()
            fq.readline()

            n_read += 1

            header_id = header[1:].split()[0]
            if "__" not in header_id:
                continue
            cb, umi = header_id.split("__", 1)

            if not seq.upper().startswith(common_prefix):
                continue
            n_prefix += 1

            variable = seq[len(common_prefix):].upper()
            tag, n_mm = best_match(variable, tag_vars, args.max_mismatch)
            if tag is None:
                continue

            out.write(f"{cb}\t{umi}\t{tag}\t{n_mm}\n")
            n_matched += 1

    print(f"reads processed:         {n_read:>10}", file=sys.stderr)
    print(f"reads with prefix:       {n_prefix:>10}", file=sys.stderr)
    print(f"reads assigned to a tag: {n_matched:>10}", file=sys.stderr)


if __name__ == "__main__":
    main()
