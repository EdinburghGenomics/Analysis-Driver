__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.util import localexecute, find_fastqs, AppLogger
import os.path

helper = TestAnalysisDriver()


def test_localexecute():
    out, err = localexecute('ls', os.path.dirname(__file__))
    print(out)
    for f in ['__init__.py', 'assets', 'test_reader', 'test_util.py', 'test_writer']:
        assert f in out.split('\n')


def test_localexecute_faulty():
    out, err = localexecute('ls', '-z', os.path.dirname(__file__))
    assert err == ('ls: illegal option -- z\n'
                   'usage: ls [-ABCFGHLOPRSTUWabcdefghiklmnopqrstuwx1] [file ...]\n')

def test_find_fastqs():
    for fq in find_fastqs(helper.fastq_path):
        assert os.path.basename(fq) in ['this.fastq.gz', 'that.fastq.gz', 'other.fastq.gz']


def test_setup_bcbio_run():
    print('Currently untestable')
    assert False


class TestLogger(TestAnalysisDriver):
    def setUp(self):
        self.logger = AppLogger()

    def test(self):
        assert True  # Not much to test here
