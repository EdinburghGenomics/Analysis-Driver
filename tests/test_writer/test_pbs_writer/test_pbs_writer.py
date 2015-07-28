__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.writer.pbs_writer.pbs_writer import PBSWriter
import os.path


class TestPBSWriter(TestAnalysisDriver):
    def setUp(self):
        self.pbs_name = 'test_' + self.__class__.__name__ + '.pbs'
        self.tmp_script = os.path.join(os.path.dirname(__file__), self.pbs_name)
        self.writer = PBSWriter(
            self.tmp_script,
            walltime='24',
            cpus='2',
            mem='1',
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    def test_init(self):
        for line in [
            '#!/bin/bash\n\n',
            '#PBS -l walltime=24:00:00\n',
            '#PBS -l ncpus=2,mem=1gb\n',
            '#PBS -N test\n',
            '#PBS -q test\n',
            '#PBS -j oe\n',
            '#PBS -o test.log\n'
        ]:
            assert line in self.writer.script

    def test_write_line(self):
        self.writer.write_line('test line')
        assert self.writer.script.endswith('test line\n')

    def test_save(self):
        self.writer.save()
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script
