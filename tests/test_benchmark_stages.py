import os
import re

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RULE_FILES = [
    os.path.join(REPO_ROOT, 'workflow', 'Snakefile'),
    os.path.join(REPO_ROOT, 'workflow', 'src', 'simulate.snmk'),
    os.path.join(REPO_ROOT, 'workflow', 'src', 'fetch_data.snmk'),
]
STAGE_TABLE = os.path.join(REPO_ROOT, 'workflow', 'data', 'benchmark_stages.tsv')

BENCHMARK_PATH = re.compile(r"'benchmarks',\s*(species \+ |f)?'([^']+)\.txt'")


def benchmark_patterns():
    """File stems of every benchmark directive, with the placeholders of the stage table."""
    patterns = set()
    for path in RULE_FILES:
        with open(path) as fh:
            for prefix, stem in BENCHMARK_PATH.findall(fh.read()):
                if prefix.startswith('species'):
                    stem = '{species}' + stem
                patterns.add(stem.replace('{_sbg_sample}', '{sample}'))
    return patterns


def stage_table_rows():
    with open(STAGE_TABLE) as fh:
        header, *rows = [line.rstrip('\n').split('\t') for line in fh if line.strip()]
    assert header == ['pattern', 'aligners', 'stage']
    return rows


def test_benchmark_directives_are_found():
    assert len(benchmark_patterns()) > 40


def test_every_benchmarked_rule_has_a_stage():
    in_table = {row[0] for row in stage_table_rows()}
    assert benchmark_patterns() - in_table == set()


def test_stage_table_has_no_stale_rule():
    in_table = {row[0] for row in stage_table_rows()}
    assert in_table - benchmark_patterns() == set()


def test_stage_table_patterns_are_unique():
    patterns = [row[0] for row in stage_table_rows()]
    assert len(patterns) == len(set(patterns))


def test_stage_table_aligners_are_known():
    known = {'starsolo', 'kallisto', 'alevin', 'sbg', 'rustody', '{aligner}'}
    for _pattern, aligners, _stage in stage_table_rows():
        assert set(aligners.split(',')) <= known
