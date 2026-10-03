from workflow.src.run_info import (
    collect_run_info,
    cpu_model,
    filesystem_type,
    total_memory_gib,
    write_run_info,
)

CPUINFO = """processor\t: 0
vendor_id\t: AuthenticAMD
model name\t: AMD EPYC 7742 64-Core Processor
cpu MHz\t\t: 2250.000

processor\t: 1
model name\t: AMD EPYC 7742 64-Core Processor
"""

MEMINFO = """MemTotal:       527990076 kB
MemFree:        12345678 kB
"""

MOUNTS = """/dev/sda1 / ext4 rw,relatime 0 0
proc /proc proc rw 0 0
/dev/sdc1 /srv/data xfs rw,noatime 0 0
server:/export /srv/data/shared nfs4 rw 0 0
"""


def test_cpu_model_reads_first_entry():
    assert cpu_model(CPUINFO) == 'AMD EPYC 7742 64-Core Processor'


def test_cpu_model_without_entry():
    assert cpu_model('processor\t: 0\nBogoMIPS\t: 50.00\n') == 'unknown'


def test_total_memory_converts_kb_to_gib():
    assert total_memory_gib(MEMINFO) == 503.5


def test_total_memory_without_entry():
    assert total_memory_gib('') == 'unknown'


def test_filesystem_type_uses_longest_mount():
    assert filesystem_type('/srv/data/shared/run1', MOUNTS) == 'nfs4'
    assert filesystem_type('/srv/data/run1', MOUNTS) == 'xfs'


def test_filesystem_type_falls_back_to_root():
    assert filesystem_type('/home/someone/output', MOUNTS) == 'ext4'


def test_filesystem_type_does_not_match_a_name_prefix():
    ## /srv/database is not under the /srv/data mount
    assert filesystem_type('/srv/database', MOUNTS) == 'ext4'


def test_filesystem_type_without_mounts():
    assert filesystem_type('/anywhere', '') == 'unknown'


def test_collect_run_info_carries_the_run_settings(tmp_path):
    info = dict(collect_run_info(str(tmp_path), cores=8, nthreads=10,
                                 max_mem_mb=20000, snakemake_version='9.0.0'))
    assert info['snakemake_cores'] == 8
    assert info['config_nthreads'] == 10
    assert info['config_max_mem_mb'] == 20000
    assert info['snakemake_version'] == '9.0.0'
    assert info['n_cpus'] >= 1
    assert set(info) >= {'cpu_model', 'memory_gib', 'kernel', 'working_dir_filesystem'}


def test_write_run_info_roundtrip(tmp_path):
    path = tmp_path / 'run_info.tsv'
    write_run_info([('cpu_model', 'some cpu'), ('n_cpus', 4)], str(path))
    assert path.read_text().splitlines() == ['key\tvalue', 'cpu_model\tsome cpu', 'n_cpus\t4']
