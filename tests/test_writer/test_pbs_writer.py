__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.writer.pbs_writer import PBSWriter
import os.path
from analysis_driver.writer.bash_commands import bcbio_env_vars, bcbio


class TestPBSWriter(TestAnalysisDriver):
    def setUp(self):
        self.pbs_name = 'test_' + self.__class__.__name__ + '.pbs'
        self.tmp_script = os.path.join(os.path.dirname(__file__), self.pbs_name)
        self.writer = PBSWriter(
            job_name='test',
            run_id='test_run',
            walltime=24,
            cpus=2,
            mem=1
        )
        self.writer.script_name = self.tmp_script

    def tearDown(self):
        try:
            os.remove(self.tmp_script)
        except FileNotFoundError:
            print('No ' + self.tmp_script + ' to remove')

    def test_init(self):
        for l in self.writer.lines:
            print(l)
        for line in [
            '#!/bin/bash\n',
            '#PBS -l walltime=24:00:00',
            '#PBS -l ncpus=2,mem=1gb',
            '#PBS -N test',
            '#PBS -q uv2000',
            '#PBS -j oe',
            '#PBS -o /Users/mwham/workspace/EdGen_Analysis_Driver/jobs/test_run/test.log'
        ]:
            assert line in self.writer.lines

    def test_write_line(self):
        self.writer.write_line('test line')
        assert self.writer.lines[-1] == 'test line'

    def test_save(self):
        print('script name')
        print(self.writer.script_name)
        print('tmp script')
        print(self.tmp_script)
        self.writer._save()
        with open(self.tmp_script) as f:
            file_lines = f.read().split('\n')
            for line in self.writer.lines:
                assert line.strip() in file_lines

    def test_bcbio(self):
        for c in bcbio_env_vars():
            self.writer.write_line(c)

        self.writer.write_line(bcbio('test.yaml', self.assets_path))

        expected = [
            '#!/bin/bash\n',
            '#PBS -l walltime=24:00:00',
            '#PBS -l ncpus=2,mem=1gb',
            '#PBS -N test',
            '#PBS -q uv2000',
            '#PBS -j oe',
            '#PBS -o /Users/mwham/workspace/EdGen_Analysis_Driver/jobs/test_run/test.log',
            '',
            'export PATH=/Users/mwham/workspace/EdGen_Analysis_Driver/Applications/tools/fake_bcbio/bin:$PATH',
            'export LD_LIBRARY_PATH=/Users/mwham/workspace/EdGen_Analysis_Driver/Applications/tools/fake_bcbio/lib:$LD_LIBRARY_PATH',
            'export PERL5LIB=/Users/mwham/workspace/EdGen_Analysis_Driver/Applications/tools/fake_bcbio/lib/perl5:$PERL5LIB',
            'export JAVA_HOME=path/to/jdk1.7.0_76/',
            'export JAVA_BINDIR=path/to/jdk1.7.0_76/bin',
            'export JAVA_ROOT=path/to/jdk1.7.0_76/',
            '',
            ('/Users/mwham/workspace/EdGen_Analysis_Driver/Applications/tools/fake_bcbio/bin/bcbio_nextgen.py'
             ' test.yaml -n 10 --workdir /Users/mwham/workspace/EdGen_Analysis_Driver/Applications/Analysis-Driver/tests/assets')
        ]
        print('\n', self.writer.lines)
        print('', expected)

        assert self.writer.lines == expected
