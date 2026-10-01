"""Machine and run settings written next to the benchmark files."""

import os
import platform


def cpu_model(cpuinfo_text):
    """First 'model name' entry of /proc/cpuinfo; some architectures have none."""
    for line in cpuinfo_text.splitlines():
        key, _, value = line.partition(':')
        if key.strip() == 'model name':
            return value.strip()
    return 'unknown'


def total_memory_gib(meminfo_text):
    """MemTotal of /proc/meminfo, which is given in kB."""
    for line in meminfo_text.splitlines():
        if line.startswith('MemTotal:'):
            return round(int(line.split()[1]) / 1024 ** 2, 1)
    return 'unknown'


def filesystem_type(path, mounts_text):
    """Type of the filesystem holding path, from the longest matching mount point."""
    path = os.path.realpath(path)
    best_mount, best_type = '', 'unknown'
    for line in mounts_text.splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        mount, fstype = fields[1], fields[2]
        inside = path == mount or path.startswith(mount.rstrip('/') + '/')
        if inside and len(mount) > len(best_mount):
            best_mount, best_type = mount, fstype
    return best_type


def _read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return ''


def collect_run_info(working_dir, cores, nthreads, max_mem_mb, snakemake_version):
    return [
        ('cpu_model', cpu_model(_read('/proc/cpuinfo'))),
        ('n_cpus', os.cpu_count()),
        ('memory_gib', total_memory_gib(_read('/proc/meminfo'))),
        ('kernel', platform.release()),
        ('working_dir_filesystem', filesystem_type(working_dir, _read('/proc/mounts'))),
        ('snakemake_version', snakemake_version),
        ('snakemake_cores', cores),
        ('config_nthreads', nthreads),
        ('config_max_mem_mb', max_mem_mb),
    ]


def write_run_info(rows, path):
    with open(path, 'w') as fh:
        fh.write('key\tvalue\n')
        for key, value in rows:
            fh.write(f'{key}\t{value}\n')
