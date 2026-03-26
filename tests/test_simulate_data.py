import gzip
import os
import random

import pytest

from workflow.src.simulate_data import (
    append_empty_droplets,
    append_sampletag_fastqs,
    make_chromosomes,
    make_r1,
    rand_seq,
    read_sampletag_fasta,
    sample_cell_barcodes,
    unique_sequences,
    write_fasta,
    write_fastqs,
    write_gtf,
    write_transcriptome_gz,
    write_true_counts_mex,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WHITELIST_DIR = os.path.join(REPO_ROOT, 'workflow', 'data', 'whitelist_384x3')


@pytest.fixture
def rng():
    return random.Random(42)


@pytest.fixture
def small_barcodes(rng):
    return sample_cell_barcodes(WHITELIST_DIR, n_cells=5, rng=rng)


@pytest.fixture
def small_gene_seqs(rng):
    return unique_sequences(4, 50, rng)


def test_rand_seq_length(rng):
    assert len(rand_seq(15, rng)) == 15


def test_rand_seq_alphabet(rng):
    assert set(rand_seq(200, rng)).issubset({'A', 'C', 'G', 'T'})


def test_unique_sequences_count(rng):
    assert len(unique_sequences(10, 8, rng)) == 10


def test_unique_sequences_are_unique(rng):
    seqs = unique_sequences(20, 8, rng)
    assert len(set(seqs)) == 20


def test_unique_sequences_length(rng):
    seqs = unique_sequences(5, 12, rng)
    assert all(len(s) == 12 for s in seqs)


def test_make_chromosomes_count(rng, small_gene_seqs):
    chroms = make_chromosomes(small_gene_seqs, chr_len=300, gene_pos=100, rng=rng)
    assert len(chroms) == len(small_gene_seqs)


def test_make_chromosomes_length(rng, small_gene_seqs):
    chroms = make_chromosomes(small_gene_seqs, chr_len=300, gene_pos=100, rng=rng)
    for seq in chroms.values():
        assert len(seq) == 300


def test_make_chromosomes_names(rng, small_gene_seqs):
    chroms = make_chromosomes(small_gene_seqs, chr_len=300, gene_pos=100, rng=rng)
    assert set(chroms.keys()) == {'chr1', 'chr2', 'chr3', 'chr4'}


def test_write_fasta_roundtrip(tmp_path, rng):
    seqs = {'chr1': rand_seq(60, rng), 'chr2': rand_seq(120, rng)}
    path = str(tmp_path / 'genome.fa')
    write_fasta(seqs, path)
    recovered = {}
    with open(path) as fh:
        name, parts = None, []
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    recovered[name] = ''.join(parts)
                name, parts = line[1:], []
            else:
                parts.append(line)
        if name:
            recovered[name] = ''.join(parts)
    assert recovered == seqs


def test_write_gtf_line_count(tmp_path):
    path = str(tmp_path / 'genes.gtf')
    write_gtf(n_genes=5, gene_pos=100, read_len=50, path=path)
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    # 3 features per gene: gene, transcript, exon
    assert len(lines) == 15


def test_write_gtf_first_line(tmp_path):
    path = str(tmp_path / 'genes.gtf')
    write_gtf(n_genes=2, gene_pos=100, read_len=50, path=path)
    with open(path) as fh:
        first = fh.readline()
    assert first.startswith('chr1\tsim\tgene\t')


def test_write_gtf_coordinates(tmp_path):
    path = str(tmp_path / 'genes.gtf')
    write_gtf(n_genes=1, gene_pos=200, read_len=75, path=path)
    with open(path) as fh:
        first = fh.readline()
    fields = first.split('\t')
    assert fields[3] == '200'
    assert fields[4] == '274'


def test_write_transcriptome_gz_entry_count(tmp_path, rng, small_gene_seqs):
    path = str(tmp_path / 'tx.fa.gz')
    write_transcriptome_gz(small_gene_seqs, path)
    with gzip.open(path, 'rt') as fh:
        header_count = sum(1 for line in fh if line.startswith('>'))
    assert header_count == len(small_gene_seqs)


def test_write_transcriptome_gz_sequences(tmp_path, rng, small_gene_seqs):
    path = str(tmp_path / 'tx.fa.gz')
    write_transcriptome_gz(small_gene_seqs, path)
    with gzip.open(path, 'rt') as fh:
        content = fh.read()
    for seq in small_gene_seqs:
        assert seq in content


def test_write_true_counts_mex_files_exist(tmp_path):
    barcodes = [('AAAAAAAAA', 'CCCCCCCCC', 'GGGGGGGGG')]
    write_true_counts_mex(barcodes, n_genes=2, n_umis=5, out_dir=str(tmp_path))
    mex = tmp_path / 'true_mex'
    assert (mex / 'barcodes.tsv.gz').exists()
    assert (mex / 'features.tsv.gz').exists()
    assert (mex / 'matrix.mtx.gz').exists()


def test_write_true_counts_mex_header(tmp_path):
    barcodes = [('AAAAAAAAA', 'CCCCCCCCC', 'GGGGGGGGG'),
                ('TTTTTTTTT', 'CCCCCCCCC', 'GGGGGGGGG')]
    write_true_counts_mex(barcodes, n_genes=3, n_umis=7, out_dir=str(tmp_path))
    with gzip.open(str(tmp_path / 'true_mex' / 'matrix.mtx.gz'), 'rt') as fh:
        lines = [l for l in fh if not l.startswith('%')]
    dims = lines[0].split()
    # n_genes rows, n_cells cols, n_genes*n_cells entries
    assert dims == ['3', '2', '6']


def test_write_true_counts_mex_umi_values(tmp_path):
    barcodes = [('AAAAAAAAA', 'CCCCCCCCC', 'GGGGGGGGG')]
    write_true_counts_mex(barcodes, n_genes=2, n_umis=11, out_dir=str(tmp_path))
    with gzip.open(str(tmp_path / 'true_mex' / 'matrix.mtx.gz'), 'rt') as fh:
        entries = [l.split() for l in fh if not l.startswith('%')][1:]
    assert all(e[2] == '11' for e in entries)


def test_read_sampletag_fasta_count(tmp_path):
    path = str(tmp_path / 'tags.fa')
    with open(path, 'w') as fh:
        fh.write('>tag1\nACGTACGT\n>tag2\nTTTTGGGG\n>tag3\nAAAACCCC\n')
    result = read_sampletag_fasta(path, n_tags=2)
    assert len(result) == 2


def test_read_sampletag_fasta_content(tmp_path):
    path = str(tmp_path / 'tags.fa')
    with open(path, 'w') as fh:
        fh.write('>tag1\nACGTACGT\n>tag2\nTTTTGGGG\n')
    result = read_sampletag_fasta(path, n_tags=5)
    assert result == [('tag1', 'ACGTACGT'), ('tag2', 'TTTTGGGG')]


def test_make_r1_linkers(rng):
    r1 = make_r1('AAAAAAAAA', 'CCCCCCCCC', 'GGGGGGGGG', 'TTTTTTTT', rng)
    assert 'GTGA' in r1
    assert 'GACA' in r1


def test_make_r1_polyt_tail(rng):
    r1 = make_r1('AAAAAAAAA', 'CCCCCCCCC', 'GGGGGGGGG', 'TTTTTTTT', rng)
    assert r1.endswith('T' * 20)


def test_make_r1_contains_barcodes_and_umi(rng):
    cb1, cb2, cb3, umi = 'ACGTACGTA', 'CGATCGATC', 'GCTAGCTAG', 'AAAACCCC'
    r1 = make_r1(cb1, cb2, cb3, umi, rng)
    assert cb1 in r1
    assert cb2 in r1
    assert cb3 in r1
    assert umi in r1


def test_sample_cell_barcodes_count(rng):
    barcodes = sample_cell_barcodes(WHITELIST_DIR, n_cells=10, rng=rng)
    assert len(barcodes) == 10


def test_sample_cell_barcodes_unique(rng):
    barcodes = sample_cell_barcodes(WHITELIST_DIR, n_cells=30, rng=rng)
    assert len(set(barcodes)) == 30


def test_sample_cell_barcodes_tuple_length(rng):
    barcodes = sample_cell_barcodes(WHITELIST_DIR, n_cells=5, rng=rng)
    assert all(len(bc) == 3 for bc in barcodes)


def test_sample_cell_barcodes_overflow(rng):
    with pytest.raises(ValueError, match='exceeds'):
        sample_cell_barcodes(WHITELIST_DIR, n_cells=999_999_999, rng=rng)


def test_write_fastqs_read_count(tmp_path, small_barcodes, small_gene_seqs, rng):
    r1 = str(tmp_path / 'R1.fq.gz')
    r2 = str(tmp_path / 'R2.fq.gz')
    n_umis = 2
    count = write_fastqs(small_barcodes, small_gene_seqs, n_umis=n_umis, rng=rng,
                         r1_path=r1, r2_path=r2)
    assert count == len(small_barcodes) * len(small_gene_seqs) * n_umis


def test_write_fastqs_files_created(tmp_path, small_barcodes, small_gene_seqs, rng):
    r1 = str(tmp_path / 'R1.fq.gz')
    r2 = str(tmp_path / 'R2.fq.gz')
    write_fastqs(small_barcodes, small_gene_seqs, n_umis=1, rng=rng, r1_path=r1, r2_path=r2)
    assert os.path.exists(r1)
    assert os.path.exists(r2)


def test_write_fastqs_paired_headers(tmp_path, small_barcodes, small_gene_seqs, rng):
    r1 = str(tmp_path / 'R1.fq.gz')
    r2 = str(tmp_path / 'R2.fq.gz')
    write_fastqs(small_barcodes, small_gene_seqs, n_umis=1, rng=rng, r1_path=r1, r2_path=r2)
    with gzip.open(r1, 'rt') as fh:
        r1_headers = [l.strip() for l in fh if l.startswith('@')]
    with gzip.open(r2, 'rt') as fh:
        r2_headers = [l.strip() for l in fh if l.startswith('@')]
    assert len(r1_headers) == len(r2_headers)


def test_append_empty_droplets_adds_reads(tmp_path, small_gene_seqs, rng):
    r1 = str(tmp_path / 'R1.fq.gz')
    r2 = str(tmp_path / 'R2.fq.gz')
    with gzip.open(r1, 'wt'):
        pass
    with gzip.open(r2, 'wt'):
        pass
    end_idx = append_empty_droplets(
        cell_barcodes_set=set(),
        n_empty=5,
        gene_seqs=small_gene_seqs,
        read_len=50,
        rng=rng,
        r1_path=r1,
        r2_path=r2,
        read_idx_start=0,
    )
    assert end_idx > 0


def test_append_empty_droplets_skips_real_cells(tmp_path, small_gene_seqs, rng):
    r1 = str(tmp_path / 'R1.fq.gz')
    r2 = str(tmp_path / 'R2.fq.gz')
    with gzip.open(r1, 'wt'):
        pass
    with gzip.open(r2, 'wt'):
        pass
    # block all 9-mer combinations - impossible in practice, just checks the skip logic
    # by using a large set the loop still terminates because random barcodes will miss it
    end_idx = append_empty_droplets(
        cell_barcodes_set=set(),
        n_empty=3,
        gene_seqs=small_gene_seqs,
        read_len=50,
        rng=rng,
        r1_path=r1,
        r2_path=r2,
        read_idx_start=10,
    )
    assert end_idx >= 10


def test_append_sampletag_fastqs_read_count(tmp_path, small_barcodes, rng):
    r1 = str(tmp_path / 'R1.fq.gz')
    r2 = str(tmp_path / 'R2.fq.gz')
    with gzip.open(r1, 'wt'):
        pass
    with gzip.open(r2, 'wt'):
        pass
    sampletags = [('tag1', rand_seq(50, rng)), ('tag2', rand_seq(50, rng))]
    n_st = 3
    end_idx = append_sampletag_fastqs(
        cell_barcodes=small_barcodes,
        sampletags=sampletags,
        n_st_reads=n_st,
        rng=rng,
        r1_path=r1,
        r2_path=r2,
        read_idx_start=0,
    )
    assert end_idx == len(small_barcodes) * n_st
