__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.writer.pbs_writer import PBSWriter
from analysis_driver.config import default as cfg
import os.path
from analysis_driver.writer.bash_commands import export_env_vars, bcbio


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
            '#PBS -o ' + cfg['jobs_dir'] + '/test_run/test.log'
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
        for c in export_env_vars():
            self.writer.write_line(c)

        self.writer.write_line(bcbio('test.yaml', self.assets_path))

        expected = [
            '#!/bin/bash\n',
            '#PBS -l walltime=24:00:00',
            '#PBS -l ncpus=2,mem=1gb',
            '#PBS -N test',
            '#PBS -q uv2000',
            '#PBS -j oe',
            '#PBS -o ' + cfg['jobs_dir'] + '/test_run/test.log',
            '#PBS -W block=true',
            '',
            'export PATH=' + cfg['bcbio'] + '/bin:$PATH',
            'export LD_LIBRARY_PATH=' + cfg['bcbio'] + '/lib:$LD_LIBRARY_PATH',
            'export PERL5LIB=' + cfg['bcbio'] + '/lib/perl5:$PERL5LIB',
            'export JAVA_HOME=' + cfg['jdk'],
            'export JAVA_BINDIR=' + cfg['jdk'] + '/bin',
            'export JAVA_ROOT=' + cfg['jdk'],
            '',
            cfg['bcbio'] + '/bin/bcbio_nextgen.py test.yaml -n 10 --workdir ' + self.assets_path
        ]
        self.compare_lists(observed=self.writer.lines, expected=expected)
