from workflow.src.cwl_step_times import main, parse_step_times, write_step_times

## cwltool 3.3 --timestamps log of a workflow with a scattered step, a nested
## workflow and a final step; colour codes as cwltool writes them to a file.
LOG = """\x1b[32m[2026-10-01 17:26:24]\x1b[0m \x1b[1;30mINFO\x1b[0m /usr/bin/cwltool 3.3.20260925135507
\x1b[32m[2026-10-01 17:26:25]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow ] start
\x1b[32m[2026-10-01 17:26:25]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow ] starting step QualCLAlign
\x1b[32m[2026-10-01 17:26:25]\x1b[0m \x1b[1;30mINFO\x1b[0m [step QualCLAlign] start
\x1b[32m[2026-10-01 17:26:25]\x1b[0m \x1b[1;30mINFO\x1b[0m [job QualCLAlign] /tmp/x$ sh \\
    -c \\
    'sleep 1; echo a' > /tmp/x/out.txt
\x1b[32m[2026-10-01 17:26:26]\x1b[0m \x1b[1;30mINFO\x1b[0m [job QualCLAlign] completed success
\x1b[32m[2026-10-01 17:26:26]\x1b[0m \x1b[1;30mINFO\x1b[0m [step QualCLAlign] start
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [job QualCLAlign_2] Max memory used: 10MiB
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [job QualCLAlign_2] completed success
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [step QualCLAlign] completed success
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow ] starting step Sub
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [step Sub] start
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow Sub] start
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow Sub] starting step Inner
\x1b[32m[2026-10-01 17:26:28]\x1b[0m \x1b[1;30mINFO\x1b[0m [step Inner] start
\x1b[32m[2026-10-01 17:26:29]\x1b[0m \x1b[1;30mINFO\x1b[0m [job Inner] completed success
\x1b[32m[2026-10-01 17:26:29]\x1b[0m \x1b[1;30mINFO\x1b[0m [step Inner] completed success
\x1b[32m[2026-10-01 17:26:29]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow Sub] completed success
\x1b[32m[2026-10-01 17:26:29]\x1b[0m \x1b[1;30mINFO\x1b[0m [step Sub] completed success
\x1b[32m[2026-10-01 17:26:29]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow ] starting step Report
\x1b[32m[2026-10-01 17:26:29]\x1b[0m \x1b[1;30mINFO\x1b[0m [step Report] start
\x1b[32m[2026-10-01 17:26:30]\x1b[0m \x1b[1;30mINFO\x1b[0m [job Report] completed success
\x1b[32m[2026-10-01 17:26:30]\x1b[0m \x1b[1;30mINFO\x1b[0m [step Report] completed success
\x1b[32m[2026-10-01 17:26:30]\x1b[0m \x1b[1;30mINFO\x1b[0m [workflow ] completed success
\x1b[32m[2026-10-01 17:26:30]\x1b[0m \x1b[1;30mINFO\x1b[0m Final process status is success
"""


def records_by_step(text=LOG):
    return {r['step']: r for r in parse_step_times(text.splitlines())}


def test_steps_are_listed_in_order_of_first_start():
    steps = [r['step'] for r in parse_step_times(LOG.splitlines())]
    assert steps == ['QualCLAlign', 'Sub', 'Inner', 'Report']


def test_scattered_step_runs_from_first_start_to_completion():
    record = records_by_step()['QualCLAlign']
    assert record['start'] == '2026-10-01 17:26:25'
    assert record['end'] == '2026-10-01 17:26:28'
    assert record['seconds'] == 3


def test_nested_workflow_step_and_its_inner_step():
    records = records_by_step()
    assert records['Sub']['seconds'] == 1
    assert records['Inner']['seconds'] == 1


def test_status_is_kept():
    assert {r['status'] for r in records_by_step().values()} == {'success'}


def test_failed_step_status():
    log = ('[2026-10-01 10:00:00] INFO [step Align] start\n'
           '[2026-10-01 10:02:05] WARNING [step Align] completed permanentFail\n')
    record = records_by_step(log)['Align']
    assert record['status'] == 'permanentFail'
    assert record['seconds'] == 125


def test_step_without_completion_has_empty_end():
    record = records_by_step('[2026-10-01 10:00:00] INFO [step Align] start\n')['Align']
    assert record['start'] == '2026-10-01 10:00:00'
    assert record['end'] == ''
    assert record['seconds'] == ''


def test_step_spanning_midnight():
    log = ('[2026-10-01 23:59:50] INFO [step Align] start\n'
           '[2026-10-02 00:00:20] INFO [step Align] completed success\n')
    assert records_by_step(log)['Align']['seconds'] == 30


def test_log_without_timestamps_gives_no_record():
    assert parse_step_times(['INFO [step Align] start', 'INFO [step Align] completed success']) == []


def test_write_step_times(tmp_path):
    path = tmp_path / 'steps.tsv'
    write_step_times(parse_step_times(LOG.splitlines()), str(path))
    lines = path.read_text().splitlines()
    assert lines[0] == 'step\tstart\tend\tseconds\tstatus'
    assert lines[1] == 'QualCLAlign\t2026-10-01 17:26:25\t2026-10-01 17:26:28\t3\tsuccess'
    assert len(lines) == 5


def test_main_reads_log_and_writes_table(tmp_path, monkeypatch):
    log = tmp_path / 'cwl.log'
    out = tmp_path / 'steps.tsv'
    log.write_text(LOG)
    monkeypatch.setattr('sys.argv', ['cwl_step_times.py', '--log', str(log), '--out', str(out)])
    main()
    assert out.read_text().splitlines()[-1].startswith('Report\t')
