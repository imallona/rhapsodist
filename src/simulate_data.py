#!/usr/bin/env python3
"""
Generates a minimal synthetic dataset for end-to-end testing of the BD Rhapsody WTA pipeline.

Outputs (all written to --out_dir):
  genome.fa                 100-chromosome genome; each chromosome contains one unique gene
                            sequence embedded at --gene_pos (1-based); flanking regions are random.
  genes.gtf                 Ensembl-style annotation; one gene / transcript / exon per chromosome.
  transcriptome.fa.gz       Gzipped FASTA of simulated transcript sequences.
  cell_barcodes.txt         Whitelist-sampled CB1+CB2+CB3 concatenations, one per simulated cell.
  sampletag_assignments.txt CB (concatenated) TAB sampletag_name, one per cell (if --sampletag_fa).
  sim_R1.fq.gz              R1 reads: BD Rhapsody dT barcode structure with random diversity inserts.
                            Headers carry ' 1:N:0:0' Illumina read-pair tag so mist_run_qualclalign.py
                            can pair files by metadata without falling back to sequence comparison.
  sim_R2.fq.gz              R2 reads: cDNA (gene sequence) + sampletag reads appended.
                            Headers carry ' 2:N:0:0' Illumina read-pair tag.

Barcode structure in R1:
  [diversity_insert (none|A|GT|TCA)] [CB1(9)] GTGA [CB2(9)] GACA [CB3(9)] [UMI(8)] [polyT(20)]

After cutadapt (min_overlap=43, noindels, e=0) and cut -c1-9,14-22,27- the standardised R1
yields CB1+CB2+CB3 (27 nt, positions 1-27) followed by UMI (8 nt, positions 28-35), which is
exactly the layout expected by the STARsolo / alevin / kallisto rules in this workflow.
"""

import argparse
import gzip
import os
import random


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--whitelist_dir', required=True,
                   help='Directory containing BD_CLS1.txt, BD_CLS2.txt, BD_CLS3.txt')
    p.add_argument('--out_dir', required=True)
    p.add_argument('--n_cells',      type=int, default=100)
    p.add_argument('--n_genes',      type=int, default=100)
    p.add_argument('--n_umis',       type=int, default=100)
    p.add_argument('--seed',         type=int, default=42)
    p.add_argument('--chr_len',      type=int, default=2000,
                   help='Total chromosome length (nt)')
    p.add_argument('--gene_pos',     type=int, default=1000,
                   help='1-based start position of the gene on each chromosome')
    p.add_argument('--read_len',     type=int, default=100,
                   help='Length of the embedded gene sequence and of cDNA R2 reads')
    p.add_argument('--sampletag_fa', default=None,
                   help='FASTA of sampletag sequences; if given, sampletag reads are added')
    p.add_argument('--n_sampletags', type=int, default=3,
                   help='Number of sampletags to use (first N from --sampletag_fa)')
    p.add_argument('--n_st_reads',   type=int, default=50,
                   help='Sampletag reads per cell')
    p.add_argument('--n_empty_droplets', type=int, default=0,
                   help='Number of empty-droplet barcodes to simulate (low-count noise for alevin knee finder)')
    return p.parse_args()


def rand_seq(length, rng):
    return ''.join(rng.choices('ACGT', k=length))


def unique_sequences(n, length, rng):
    seen, seqs = set(), []
    while len(seqs) < n:
        s = rand_seq(length, rng)
        if s not in seen:
            seen.add(s)
            seqs.append(s)
    return seqs


def make_chromosomes(gene_seqs, chr_len, gene_pos, rng):
    gene_len = len(gene_seqs[0])
    chroms = {}
    for idx, gseq in enumerate(gene_seqs):
        pre  = rand_seq(gene_pos - 1, rng)
        post = rand_seq(chr_len - (gene_pos - 1) - gene_len, rng)
        chroms[f'chr{idx + 1}'] = pre + gseq + post
    return chroms


def write_fasta(seqs_dict, path, line_width=60):
    with open(path, 'w') as fh:
        for name, seq in seqs_dict.items():
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + '\n')


def write_gtf(n_genes, gene_pos, read_len, path):
    width = len(str(n_genes))
    with open(path, 'w') as fh:
        for i in range(n_genes):
            chrom = f'chr{i + 1}'
            gid   = f'gene{i + 1:0{width}d}'
            tid   = f'{gid}_tx'
            start = gene_pos
            end   = gene_pos + read_len - 1
            gene_attr = f'gene_id "{gid}"; gene_version "1";'
            tx_attr   = f'{gene_attr} transcript_id "{tid}"; transcript_version "1";'
            exon_attr = (f'{tx_attr} exon_number "1";'
                         f' exon_id "{gid}_exon1"; exon_version "1";')
            for feat, attr in (('gene', gene_attr),
                                ('transcript', tx_attr),
                                ('exon', exon_attr)):
                fh.write(f'{chrom}\tsim\t{feat}\t{start}\t{end}\t.\t+\t.\t{attr}\n')


def write_transcriptome_gz(gene_seqs, path):
    width = len(str(len(gene_seqs)))
    with gzip.open(path, 'wt') as fh:
        for i, seq in enumerate(gene_seqs):
            gid = f'gene{i + 1:0{width}d}'
            tid = f'{gid}_tx'
            fh.write(f'>{tid} {gid}\n{seq}\n')


def write_true_counts_mex(cell_barcodes, n_genes, n_umis, out_dir):
    """Write ground-truth count matrix (MEX format) to out_dir/true_mex/."""
    import struct
    mex_dir = os.path.join(out_dir, 'true_mex')
    os.makedirs(mex_dir, exist_ok=True)
    n_cells = len(cell_barcodes)
    width = len(str(n_genes))

    with gzip.open(os.path.join(mex_dir, 'barcodes.tsv.gz'), 'wt') as fh:
        for cb1, cb2, cb3 in cell_barcodes:
            fh.write(f'{cb1}{cb2}{cb3}\n')

    with gzip.open(os.path.join(mex_dir, 'features.tsv.gz'), 'wt') as fh:
        for i in range(n_genes):
            gid = f'gene{i + 1:0{width}d}'
            fh.write(f'{gid}\t{gid}\tGene Expression\n')

    n_entries = n_cells * n_genes
    with gzip.open(os.path.join(mex_dir, 'matrix.mtx.gz'), 'wt') as fh:
        fh.write('%%MatrixMarket matrix coordinate integer general\n%\n')
        fh.write(f'{n_genes} {n_cells} {n_entries}\n')
        for j in range(1, n_cells + 1):
            for i in range(1, n_genes + 1):
                fh.write(f'{i} {j} {n_umis}\n')


def sample_cell_barcodes(whitelist_dir, n_cells, rng):
    wl = []
    for cls in ('BD_CLS1.txt', 'BD_CLS2.txt', 'BD_CLS3.txt'):
        with open(os.path.join(whitelist_dir, cls)) as fh:
            wl.append([line.strip() for line in fh if line.strip()])
    max_cells = len(wl[0]) * len(wl[1]) * len(wl[2])
    if n_cells > max_cells:
        raise ValueError(
            f"Requested n_cells={n_cells} exceeds the {max_cells} unique barcode "
            "combinations available from the provided whitelists."
        )
    chosen, barcodes = set(), []
    while len(barcodes) < n_cells:
        combo = (rng.choice(wl[0]), rng.choice(wl[1]), rng.choice(wl[2]))
        if combo not in chosen:
            chosen.add(combo)
            barcodes.append(combo)
    return barcodes


def read_sampletag_fasta(path, n_tags):
    tags, name, seq_parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    tags[name] = ''.join(seq_parts)
                name = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            tags[name] = ''.join(seq_parts)
    names = list(tags)[:n_tags]
    return [(n, tags[n]) for n in names]


def make_r1(cb1, cb2, cb3, umi, rng):
    diversity_inserts = ['', 'A', 'GT', 'TCA']
    div = rng.choice(diversity_inserts)
    return div + cb1 + 'GTGA' + cb2 + 'GACA' + cb3 + umi + 'T' * 20


def write_fastqs(cell_barcodes, gene_seqs, n_umis, rng, r1_path, r2_path):
    read_idx = 0
    with gzip.open(r1_path, 'wt') as fq1, gzip.open(r2_path, 'wt') as fq2:
        for cb1, cb2, cb3 in cell_barcodes:
            for gseq in gene_seqs:
                umis = unique_sequences(n_umis, 8, rng)
                for umi in umis:
                    r1 = make_r1(cb1, cb2, cb3, umi, rng)
                    tag = f'sim{read_idx}'
                    fq1.write(f'@{tag} 1:N:0:0\n{r1}\n+\n{"I" * len(r1)}\n')
                    fq2.write(f'@{tag} 2:N:0:0\n{gseq}\n+\n{"I" * len(gseq)}\n')
                    read_idx += 1
    return read_idx


def append_empty_droplets(cell_barcodes_set, n_empty, gene_seqs, read_len, rng,
                          r1_path, r2_path, read_idx_start):
    """Append ambient-RNA reads from barcodes NOT in the true cell set.

    Each empty droplet gets 1-10 reads from a random gene, creating the
    bimodal barcode-count distribution that alevin's knee finder requires.
    """
    read_idx = read_idx_start
    added = 0
    with gzip.open(r1_path, 'at') as fq1, gzip.open(r2_path, 'at') as fq2:
        while added < n_empty:
            # random barcode tuple unlikely to collide with real cells
            cb1 = rand_seq(9, rng)
            cb2 = rand_seq(9, rng)
            cb3 = rand_seq(9, rng)
            if (cb1, cb2, cb3) in cell_barcodes_set:
                continue
            n_reads = rng.randint(1, 10)
            gseq = rng.choice(gene_seqs)
            for _ in range(n_reads):
                umi = rand_seq(8, rng)
                r1 = make_r1(cb1, cb2, cb3, umi, rng)
                tag = f'emp{read_idx}'
                fq1.write(f'@{tag} 1:N:0:0\n{r1}\n+\n{"I" * len(r1)}\n')
                fq2.write(f'@{tag} 2:N:0:0\n{gseq}\n+\n{"I" * len(gseq)}\n')
                read_idx += 1
            added += 1
    return read_idx


def append_sampletag_fastqs(cell_barcodes, sampletags, n_st_reads, rng,
                             r1_path, r2_path, read_idx_start):
    read_idx = read_idx_start
    n_tags = len(sampletags)
    with gzip.open(r1_path, 'at') as fq1, gzip.open(r2_path, 'at') as fq2:
        for cell_idx, (cb1, cb2, cb3) in enumerate(cell_barcodes):
            _st_name, st_seq = sampletags[cell_idx % n_tags]
            umis = unique_sequences(n_st_reads, 8, rng)
            for umi in umis:
                r1 = make_r1(cb1, cb2, cb3, umi, rng)
                tag = f'st{read_idx}'
                fq1.write(f'@{tag} 1:N:0:0\n{r1}\n+\n{"I" * len(r1)}\n')
                fq2.write(f'@{tag} 2:N:0:0\n{st_seq}\n+\n{"I" * len(st_seq)}\n')
                read_idx += 1
    return read_idx


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = random.Random(args.seed)

    gene_seqs = unique_sequences(args.n_genes, args.read_len, rng)
    chroms    = make_chromosomes(gene_seqs, args.chr_len, args.gene_pos, rng)

    write_fasta(chroms,   os.path.join(args.out_dir, 'genome.fa'))
    write_gtf(args.n_genes, args.gene_pos, args.read_len,
              os.path.join(args.out_dir, 'genes.gtf'))
    write_transcriptome_gz(gene_seqs, os.path.join(args.out_dir, 'transcriptome.fa.gz'))

    cell_barcodes = sample_cell_barcodes(args.whitelist_dir, args.n_cells, rng)
    with open(os.path.join(args.out_dir, 'cell_barcodes.txt'), 'w') as fh:
        for cb1, cb2, cb3 in cell_barcodes:
            fh.write(f'{cb1}{cb2}{cb3}\n')

    next_idx = write_fastqs(cell_barcodes, gene_seqs, args.n_umis, rng,
                            os.path.join(args.out_dir, 'sim_R1.fq.gz'),
                            os.path.join(args.out_dir, 'sim_R2.fq.gz'))

    write_true_counts_mex(cell_barcodes, args.n_genes, args.n_umis, args.out_dir)

    print(f'Generated {args.n_genes} chromosomes, {args.n_cells} cells x '
          f'{args.n_genes} genes x {args.n_umis} UMIs = {next_idx} cDNA reads')

    if args.n_empty_droplets > 0:
        cell_barcodes_set = set(cell_barcodes)
        next_idx = append_empty_droplets(cell_barcodes_set, args.n_empty_droplets,
                                         gene_seqs, args.read_len, rng,
                                         os.path.join(args.out_dir, 'sim_R1.fq.gz'),
                                         os.path.join(args.out_dir, 'sim_R2.fq.gz'),
                                         next_idx)
        print(f'Added {args.n_empty_droplets} empty-droplet barcodes')

    if args.sampletag_fa:
        sampletags = read_sampletag_fasta(args.sampletag_fa, args.n_sampletags)
        append_sampletag_fastqs(cell_barcodes, sampletags, args.n_st_reads, rng,
                                os.path.join(args.out_dir, 'sim_R1.fq.gz'),
                                os.path.join(args.out_dir, 'sim_R2.fq.gz'),
                                next_idx)
        with open(os.path.join(args.out_dir, 'sampletag_assignments.txt'), 'w') as fh:
            fh.write('cell_barcode\tsampletag\n')
            for cell_idx, (cb1, cb2, cb3) in enumerate(cell_barcodes):
                st_name, _ = sampletags[cell_idx % len(sampletags)]
                fh.write(f'{cb1}{cb2}{cb3}\t{st_name}\n')
        print(f'Added {args.n_cells * args.n_st_reads} sampletag reads '
              f'({args.n_st_reads}/cell across {len(sampletags)} tags)')


if __name__ == '__main__':
    main()
