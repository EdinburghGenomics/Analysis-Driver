__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.writer import bash_commands
from analysis_driver.reader import SampleSheet, RunInfo
from analysis_driver.config import default as cfg
import os.path

helper = TestAnalysisDriver()
sample_sheet_csv = os.path.join(helper.assets_path, 'SampleSheet_analysis_driver.csv')
sample_sheet = SampleSheet(helper.assets_path)
run_info = RunInfo(helper.assets_path)


def test_bcl2fastq():
    mask = sample_sheet.generate_mask(run_info.mask)
    expected = (
        cfg['bcl2fastq'] +
        ' -l INFO'
        ' --runfolder-dir ' + helper.assets_path +
        ' --output-dir ' + helper.fastq_path +
        ' -r 8 -d 8 -p 8 -w 8' +
        ' --sample-sheet ' + sample_sheet_csv +
        ' --use-bases-mask ' + mask
    )
    observed = bash_commands.bcl2fastq(helper.assets_path, helper.fastq_path, sample_sheet_csv, mask)

    assert observed == expected


def test_fastqc():
    test_fastq = os.path.join(helper.fastq_path, '10015AT', '10015ATA0001L05', 'this.fastq.gz')
    expected = '--nogroup -t 4 -q ' + test_fastq
    assert bash_commands.fastqc(test_fastq).endswith(expected)


def test_bcbio():
    run_yaml = os.path.join(helper.assets_path, 'run.yaml')
    expected = '%s -n 16 --workdir %s' % (run_yaml, helper.assets_path)
    assert bash_commands.bcbio(run_yaml, helper.assets_path, threads=16).endswith(expected)
