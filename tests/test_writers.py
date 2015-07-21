__author__ = 'mwham'

import os

from tests import TestAnalysisDriver

from analysis_driver.writer import BCBioCSVWriter
from analysis_driver.reader.sample_sheet import SampleSheet

from analysis_driver.writer.pbs_writer.pbs_writer import PBSWriter
from analysis_driver.writer.pbs_writer.bcbio_pbs_writer import BCBioPBSWriter
from analysis_driver.writer.pbs_writer.bcl2fastq_pbs_writer import BCL2FastqPBSWriter
from analysis_driver.writer.pbs_writer.fastqc_pbs_writer import FastqcPBSWriter

# TODO: move these back to separate files


class TestPBSWriter(TestAnalysisDriver):
    def __init__(self):
        super().__init__()
        self.tmp_script = os.path.join(os.path.dirname(__file__), self.pbs_name)

    def setUp(self):
        self.writer = PBSWriter(
            self.tmp_script,
            walltime='24',
            cpus='2',
            mem='1',
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    @property
    def pbs_name(self):
        return 'test_' + self.__class__.__name__ + '.pbs'

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

    def tearDown(self):
        os.remove(self.tmp_script)


class TestBCBioCSVWriter(TestAnalysisDriver):
    def setUp(self):
        sample_sheet = SampleSheet(self.assets_path)
        self.tmp_run_dir = os.path.join(os.path.dirname(__file__), 'test_run_dir')
        self.writer = BCBioCSVWriter(self.fastq_path, self.tmp_run_dir, sample_sheet)

    def test_write(self):
        self.writer.write()
        with open(self.tmp_run_dir) as f:
            assert True
            # TODO: assert that the file has correct content

    def test_find_fastqs(self):
        fastqs = self.writer._find_fastqs(self.fastq_path, '10015AT')
        for file_name in ['this.fastq', 'that.fastq', 'other.fastq']:
            assert os.path.join(self.fastq_path, '10015AT', file_name) in fastqs


class TestBCBioPBSWriter(TestPBSWriter):
    bcbio_path = os.path.join('/', 'path', 'to', 'bcbio_nextgen.py')
    run_yaml = os.path.join('/', 'path', 'to', 'bcbio_run.yaml')
    workdir = os.path.join('/', 'path', 'to', 'bcbio', 'work')

    def __init__(self):
        super().__init__()

    def setUp(self):
        self.writer = BCBioPBSWriter('test', 'test_bcbio', 'test.log')

    def test_init(self):
        for line in [
            '#!/bin/bash\n\n',
            '#PBS -l walltime=72:00:00\n',
            '#PBS -l ncpus=8,mem=64gb\n',
            '#PBS -N test\n',
            '#PBS -q test\n',
            '#PBS -j oe\n',
            '#PBS -o test.log\n'
        ]:
            assert line in self.writer.script

    def test_bcbio(self):
        self.writer._bcbio(
            self.bcbio_path,
            self.run_yaml,
            self.workdir
        )
        for line in [
            'export JAVA_HOME=/home/U008/edingen/Applications/jdk1.7.0_76/',
            'export JAVA_BINDIR=/home/U008/edingen/Applications/jdk1.7.0_76/bin',
            'export JAVA_ROOT=/home/edingen/Applications/jdk1.7.0_76/\n',
            '%s %s -n 16 --workdir %s\n' % (self.bcbio_path, self.run_yaml, self.workdir)
        ]:
            assert line in self.writer.script

    def test_write(self):
        self.writer.write(self.bcbio_path, self.run_yaml, self.workdir)
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script


class TestBCL2FastqPBSWriter(TestPBSWriter):
    def __init__(self):
        super().__init__()

    def setUp(self):
        self.writer = BCL2FastqPBSWriter(
            self.tmp_script,
            walltime='24',
            cpus='2',
            mem='1',
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    def test_bcl2fastq(self):
        self.writer._bcl2fastq(mask='this,that,other', input_dir=self.assets_path, fastq_path='a_fastq_path')
        assert self.writer.script.endswith(
            'bcl2fastq -l INFO --runfolder-dir %s --output-dir %s --sample-sheet %s --use-bases-mask %s\n' % (
                self.assets_path,
                'a_fastq_path',
                os.path.join(self.assets_path, 'SampleSheet.csv'),
                'this,that,other'
            )
        )

    def test_write(self):
        self.writer.write(mask='this,that,other', input_dir=self.assets_path, fastq_path='a_fastq_path')
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script


class TestFastqcPBSWriter(TestPBSWriter):
    def __init__(self):
        super().__init__()

    def setUp(self):
        self.writer = FastqcPBSWriter(
            self.tmp_script,
            walltime='24',
            cpus='2',
            mem='1',
            job_name='test',
            log_file='test.log',
            queue='test'
        )

    def test_fastqc(self):
        self.writer._fastqc(self.assets_path)
        for line in ['FASTQ_FILES=`find ' + self.assets_path + ' -name \'*.fastq.gz\'`\n',
                     'fastqc --nogroup -t 8 -q $FASTQ_FILES\n']:
            assert line in self.writer.script

    def test_write(self):
        self.writer.write(input_dir=self.assets_path, run_dir=self.assets_path)
        with open(self.tmp_script) as f:
            assert f.read() == self.writer.script
