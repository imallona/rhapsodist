#!/usr/bin/env python3
##
## cli wrapper around the rhapsodist snakemake pipeline.
##
## usage:
##   rhapsodist --configfile config.yaml [--cores n] [extra snakemake args ...]
##
## all unrecognised arguments are forwarded to snakemake as-is.
##
## izaskun mallona

import argparse
import os
import subprocess
import sys
from pathlib import Path


def find_snakefile():
    """locate snakefile: next to this file's package root, then cwd."""
    candidate = Path(__file__).resolve().parent.parent / "Snakefile"
    if candidate.exists():
        return candidate
    candidate = Path("Snakefile")
    if candidate.exists():
        return candidate
    sys.exit(
        "cannot find Snakefile. run rhapsodist from the repository root "
        "or install the package with 'pip install -e .'"
    )


def main():
    ap = argparse.ArgumentParser(
        prog="rhapsodist",
        description="run the rhapsodist bd rhapsody wta pipeline.",
        epilog="any extra arguments are forwarded verbatim to snakemake.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--configfile",
        default="configs/config.yaml",
        metavar="FILE",
        help="path to the pipeline config yaml (default: configs/config.yaml)",
    )
    ap.add_argument(
        "--cores",
        default=str(os.cpu_count() or 1),
        metavar="N",
        help="number of cores to pass to snakemake (default: all available)",
    )
    ap.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="pass --dry-run to snakemake",
    )

    args, extra = ap.parse_known_args()

    snakefile = find_snakefile()

    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", args.configfile,
        "--cores", args.cores,
        "--use-conda",
    ]
    if args.dry_run:
        cmd.append("--dry-run")
    cmd.extend(extra)

    print("running:", " ".join(cmd), file=sys.stderr)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
