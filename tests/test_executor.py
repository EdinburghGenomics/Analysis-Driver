import pytest
import os
import shutil
import sys
import logging
import subprocess
from unittest.mock import patch, Mock
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.executor import script_writers, Executor, StreamExecutor, ArrayExecutor
from analysis_driver.executor.executor import ClusterExecutor, PBSExecutor, SlurmExecutor
from analysis_driver.exceptions import AnalysisDriverError

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

get_stdout = 'analysis_driver.executor.executor.ClusterExecutor._get_stdout'
sleep = 'analysis_driver.executor.executor.sleep'


class TestScriptWriter(TestAnalysisDriver):
    array_index = 'JOB_INDEX'
    exp_header = []

    def setUp(self):
        self.working_dir = os.path.join(self.assets_path, 'test_script_writer_wd')
        os.makedirs(self.working_dir, exist_ok=True)
        self.script_writer = script_writers.ScriptWriter('a_job_name', self.working_dir, 'a_job_queue')

    def tearDown(self):
        shutil.rmtree(self.working_dir)

    def _compare_writer_lines(self, expected):
        self.compare_lists(
            [l.rstrip('\n') for l in self.script_writer.lines],
            [l.rstrip('\n') for l in self.exp_header + expected]
        )

    def test_write_line(self):
        self.script_writer.write_line('a_line')
        self._compare_writer_lines(['a_line'])

    def test_write_job(self):
        self.script_writer.write_jobs(['a_cmd'])
        self._compare_writer_lines(['a_cmd'])

    def test_write_job_prelim_cmds(self):
        self.script_writer.write_jobs(['a_cmd'], prelim_cmds=['a_prelim_cmd'])
        self._compare_writer_lines(['a_prelim_cmd', '', 'a_cmd'])

    def test_start_array(self):
        self.script_writer._start_array()
        self._compare_writer_lines(['case ${array_index} in'.format(array_index=self.array_index)])

    def test_finish_array(self):
        self.script_writer._finish_array()
        self._compare_writer_lines(
            [
                '*) echo "Unexpected {array_index}: ${array_index}"'.format(array_index=self.array_index),
                'esac'
            ]
        )

    def test_write_array_cmd(self):
        self.script_writer._write_array_cmd(1337, 'an_array_cmd')
        self.script_writer._write_array_cmd(
            1338, 'another_array_cmd', log_file=os.path.join(self.assets_path, 'a_log_file')
        )
        self._compare_writer_lines(
            [
                '1337) an_array_cmd\n' + ';;',
                '1338) another_array_cmd > ' + os.path.join(self.assets_path, 'a_log_file') + ' 2>&1''\n;;'
            ]
        )

    def test_write_job_array(self):
        self.script_writer.write_jobs(['a_cmd', 'another_cmd'])
        expected = [
            'case ${array_index} in'.format(array_index=self.array_index),
            '1) a_cmd > ' + self.script_writer.log_file + '1 2>&1' + '\n' + ';;',
            '2) another_cmd > ' + self.script_writer.log_file + '2 2>&1' + '\n' + ';;',
            '*) echo "Unexpected {array_index}: ${array_index}"'.format(array_index=self.array_index),
            'esac'
        ]
        self._compare_writer_lines(expected)

    def test_write_job_array_prelim_cmds(self):
        self.script_writer.write_jobs(['a_cmd', 'another_cmd'], prelim_cmds=['a_prelim_cmd'])
        expected = [
            'a_prelim_cmd',
            '',
            'case ${array_index} in'.format(array_index=self.array_index),
            '1) a_cmd > ' + self.script_writer.log_file + '1 2>&1' + '\n' + ';;',
            '2) another_cmd > ' + self.script_writer.log_file + '2 2>&1' + '\n' + ';;',
            '*) echo "Unexpected {array_index}: ${array_index}"'.format(array_index=self.array_index),
            'esac'
        ]
        self._compare_writer_lines(expected)

    def test_save(self):
        self.script_writer.write_line('a_line')
        self.script_writer._save()
        assert open(self.script_writer.script_name, 'r').readlines() == ['a_line\n']

    def test_trim_field(self):
        assert self.script_writer._trim_field('a_field_name_too_long_for_pbs', 15) == 'a_field_name_to'


class TestPBSWriter(TestScriptWriter):
    array_index = 'PBS_ARRAY_INDEX'

    def setUp(self):
        super().setUp()
        self.script_writer = script_writers.PBSWriter(
            'a_job_name',
            self.working_dir,
            'a_job_queue',
            walltime=3,
            cpus=2,
            mem=1,
            jobs=1,
            log_commands=True
        )
        self.exp_header = [
            '#!/bin/bash\n',
            '#PBS -l ncpus=2,mem=1gb',
            '#PBS -q a_job_queue',
            '#PBS -j oe',
            '#PBS -o ' + os.path.join(self.working_dir, 'a_job_name.log'),
            '#PBS -W block=true',
            '#PBS -l walltime=3:00:00',
            '#PBS -N a_job_name',
            'cd ' + self.script_writer.working_dir,
            ''
        ]

    def test_write_header(self):
        self.compare_lists(self.script_writer.lines, self.exp_header)

    def test_write_header_no_walltime(self):
        script_writer = script_writers.PBSWriter(
            'a_job_name',
            self.working_dir,
            'a_job_queue',
            walltime=None,
            cpus=2,
            mem=1,
            jobs=1,
            log_commands=True
        )
        exp_header = [
            '#!/bin/bash\n',
            '#PBS -l ncpus=2,mem=1gb',
            '#PBS -q a_job_queue',
            '#PBS -j oe',
            '#PBS -o ' + os.path.join(self.working_dir, 'a_job_name.log'),
            '#PBS -W block=true',
            '#PBS -N a_job_name',
            'cd ' + self.script_writer.working_dir,
            ''
        ]
        self.compare_lists(script_writer.lines, exp_header)

    def test_start_array(self):
        self.script_writer._start_array()
        self._compare_writer_lines(['case $PBS_ARRAY_INDEX in\n'])

    def test_save(self):
        self.script_writer.write_line('a_line')
        self.script_writer._save()
        expected = '\n'.join(self.exp_header + ['a_line', ''])
        assert open(self.script_writer.script_name, 'r').read() == expected


class TestExecutor(TestAnalysisDriver):
    def _get_executor(self, cmd):
        return Executor(cmd)

    def test_cmd(self):
        e = self._get_executor('ls ' + os.path.join(self.assets_path, '..'))
        exit_status = e.join()
        assert exit_status == 0

    def test_dodgy_cmd(self):
        with pytest.raises(AnalysisDriverError) as err:
            e = self._get_executor('dodgy_cmd')
            e.join()
            assert 'Command failed: \'dodgy_cmd\'' in str(err)

    def test_process(self):
        e = self._get_executor('ls ' + os.path.join(self.assets_path, '..'))
        assert e.proc is None
        proc = e._process()
        assert proc is e.proc and isinstance(e.proc, subprocess.Popen)


class TestStreamExecutor(TestExecutor):
    def _get_executor(self, cmd):
        return StreamExecutor(cmd)

    def test_cmd(self):
        e = self._get_executor(os.path.join(self.assets_path, 'countdown.sh'))
        e.start()
        assert e.join() == 0

    def test_dodgy_command(self):
        e = self._get_executor(os.path.join(self.assets_path, 'countdown.sh') + ' dodgy')
        e.start()
        assert e.join() == 13  # same exit status as the running script

    def test_dodgy_cmd(self):
        with pytest.raises(AnalysisDriverError) as err:
            e = self._get_executor('dodgy_cmd')
            e.start()
            e.join()
            assert 'self.proc command failed: \'dodgy_cmd\'' in str(err)


class TestArrayExecutor(TestExecutor):
    def _get_executor(self, cmds):
        return ArrayExecutor(cmds, stream=True)

    def test_cmd(self):
        e = self._get_executor(['ls', 'ls -lh', 'pwd'])
        e.start()
        assert e.join() == 0
        assert e.exit_statuses == [0, 0, 0]

    def test_dodgy_cmd(self):
        e = self._get_executor(['ls', 'non_existent_cmd', 'pwd'])
        e.start()
        with pytest.raises(AnalysisDriverError) as err:
            e.join()
            assert 'Commands failed' in str(err)


class TestClusterExecutor(TestAnalysisDriver):
    e_cls = ClusterExecutor

    @property
    def script(self):
        return os.path.join(self.assets_path, 'countdown.sh')

    def setUp(self):
        os.makedirs(os.path.join(self.assets_path, 'a_run_id'), exist_ok=True)
        self.executor = self._get_executor(self.script)

    def tearDown(self):
        shutil.rmtree(os.path.join(self.assets_path, 'a_run_id'))

    def _get_executor(self, cmd):
        get_writer = 'analysis_driver.executor.executor.ClusterExecutor._get_writer'
        with patch(get_writer, return_value=Mock(script_name='a_script_name')):
            return self.e_cls(
                cmd,
                job_name='test_job',
                working_dir=os.path.join(self.assets_path, 'a_run_id')
            )

    def test_get_stdout(self):
        popen = 'analysis_driver.executor.executor.subprocess.Popen'
        with patch(popen, return_value=Mock(wait=Mock(return_value=None))) as p:
            assert self.executor._get_stdout('ls -d ' + self.assets_path).endswith('tests/assets')
            p.assert_called_with(['ls', '-d', self.assets_path], stdout=-1, stderr=-1)

    def test_cmd(self):
        assert self.executor.cmd == '/bin/sh a_script_name'

    @patch(get_stdout, return_value=None)
    def test_dodgy_cmd(self, mocked_get_stdout):
        with pytest.raises(AnalysisDriverError) as e:
            ex = self._get_executor(os.path.join(self.assets_path, 'non_existent_script.sh'))
            ex.cmd = '/bin/sh non_existent_script.sh'
            ex.start()
            assert str(e) == 'Job submission failed'
        mocked_get_stdout.assert_called_with('/bin/sh non_existent_script.sh')

    def test_join(self):
        job_finished = 'analysis_driver.executor.executor.' + self.e_cls.__name__ + '._job_finished'
        exit_code = 'analysis_driver.executor.executor.' + self.e_cls.__name__ + '._job_exit_code'
        self.executor.finished_statuses = 'FXM'
        with patch(job_finished, return_value=True), patch(exit_code, return_value=0), patch(sleep):
            assert self.executor.join() == 0


class TestPBSExecutor(TestClusterExecutor):
    e_cls = PBSExecutor

    def test_qstat(self):
        with patch(get_stdout, return_value='this\nthat\nother') as p:
            assert self.executor._qstat() == 'other'.split()
            p.assert_called_with('qstat -x None')

    def test_job_status(self):
        qstat = 'analysis_driver.executor.executor.PBSExecutor._qstat'
        fake_report = ('1337', 'a_job', 'a_user', '10:00:00', 'R',  'q')
        with patch(qstat, return_value=fake_report) as p:
            assert self.executor._job_status() == 'R'

    def test_job_finished(self):
        job_status = 'analysis_driver.executor.executor.PBSExecutor._qstat'
        with patch(job_status, return_value=(1, '1', 'user', 'time', 'B', 'queue')):
            assert not self.executor._job_finished()
        with patch(job_status, return_value=(1, '1', 'user', 'time', 'F', 'queue')):
            assert self.executor._job_finished()


class TestSlurmExecutor(TestClusterExecutor):
    e_cls = SlurmExecutor

    def test_sacct(self):
        with patch(get_stdout, return_value='1:0') as p:
            assert self.executor._sacct('ExitCode') == '1:0'
            p.assert_called_with('sacct -n -j None -o ExitCode')

    def test_job_finished(self):
        sacct = 'analysis_driver.executor.executor.SlurmExecutor._sacct'
        patched_squeue = patch('analysis_driver.executor.executor.SlurmExecutor._squeue', return_value='')
        with patch(sacct, return_value='RUNNING'), patched_squeue:
            assert not self.executor._job_finished()
        with patch(sacct, return_value='COMPLETED'), patched_squeue:
            assert self.executor._job_finished()

    def test_job_exit_code(self):
        sacct = 'analysis_driver.executor.executor.SlurmExecutor._sacct'
        with patch(sacct, return_value='CANCELLED 0:0'):
            assert self.executor._job_exit_code() == 9
        with patch(sacct, return_value='COMPLETED 0:x'):
            assert self.executor._job_exit_code() == 0
