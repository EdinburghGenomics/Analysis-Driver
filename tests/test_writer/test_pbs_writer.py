__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.writer.pbs_writer import PBSWriter
import os.path


class TestPBSWriter(TestAnalysisDriver):
    def setUp(self):
        self.pbs_name = 'test_' + self.__class__.__name__ + '.pbs'
        self.tmp_script = os.path.join(os.path.dirname(__file__), self.pbs_name)
        self.writer = PBSWriter(
            self.tmp_script,
            walltime=24,
            cpus=2,
            mem=1,
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    def tearDown(self):
        try:
            os.remove(self.tmp_script)
        except FileNotFoundError:
            print('No ' + self.tmp_script + ' to remove')

    def test_init(self):
        for line in [
            '#!/bin/bash\n',
            '#PBS -l walltime=24:00:00',
            '#PBS -l ncpus=2,mem=1gb',
            '#PBS -N test',
            '#PBS -q test',
            '#PBS -j oe',
            '#PBS -o test.log'
        ]:
            assert line in self.writer.lines

    def test_write_line(self):
        self.writer.write_line('test line')
        assert self.writer.lines[-1] == 'test line'

    def test_save(self):
        self.writer.save()
        with open(self.tmp_script) as f:
            file_lines = f.read().split('\n')
            for line in self.writer.lines:
                assert line.strip() in file_lines
