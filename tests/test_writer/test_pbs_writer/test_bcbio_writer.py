__author__ = 'mwham'
from tests.test_writer.test_pbs_writer.test_pbs_writer import TestPBSWriter
from analysis_driver.writer.pbs_writer.bcbio_writer import BCBioWriter
import os.path


class TestBCBioWriter(TestPBSWriter):
    bcbio_path = os.path.join('/', 'home', 'U008', 'edingen', 'Applications', 'bcbio', 'bin', 'bcbio_nextgen.py')
    run_yaml = os.path.join('/', 'path', 'to', 'bcbio_run.yaml')
    workdir = os.path.join('/', 'path', 'to', 'bcbio', 'work')

    def setUp(self):
        super().setUp()
        self.writer = BCBioWriter(self.tmp_script, 'test_bcbio', 'test.log', 1, queue='test')

    def test_init(self):
        for line in [
            '#!/bin/bash\n\n',
            '#PBS -l walltime=72:00:00\n',
            '#PBS -l ncpus=8,mem=64gb\n',
            '#PBS -N test_bcbio\n',
            '#PBS -q test\n',
            '#PBS -j oe\n',
            '#PBS -o test.log\n'
        ]:
            assert line in self.writer.script

    def test_bcbio(self):
        self.writer.add_bcbio_job(
            self.run_yaml,
            self.workdir
        )
        print(self.writer.script)
        for line in [
            'export JAVA_HOME=/home/U008/edingen/Applications/jdk1.7.0_76/',
            'export JAVA_BINDIR=/home/U008/edingen/Applications/jdk1.7.0_76/bin',
            'export JAVA_ROOT=/home/edingen/Applications/jdk1.7.0_76/\n',
            '%s %s -n 16 --workdir %s\n' % (self.bcbio_path, self.run_yaml, self.workdir)
        ]:

            assert line in self.writer.script

    def test_write(self):
        self.writer.add_bcbio_job(self.run_yaml, self.workdir)
        self.writer.write()
        self.writer.pbs_file.close()
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script

    def test_save(self):
        pass
