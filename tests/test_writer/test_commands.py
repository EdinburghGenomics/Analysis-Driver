__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.writer import commands
from analysis_driver.reader import SampleSheet
import os.path

helper = TestAnalysisDriver()
sample_sheet = SampleSheet(helper.assets_path)


def test_bcl2fastq():
    input_dir = helper.assets_path
    fastq_path = helper.fastq_path
    mask = sample_sheet.generate_mask()
    expected = ('-l INFO --runfolder-dir %s '
                '--output-dir %s --sample-sheet %s '
                '--use-bases-mask %s '
                '-r 8 -d 8 -p 8 -w 8') % (
        input_dir,
        fastq_path,
        os.path.join(input_dir, 'SampleSheet_analysis_driver.csv'),
        mask
    )

    assert commands.bcl2fastq(mask, input_dir, fastq_path).endswith(expected)


def test_fastqc():
    test_fastq = os.path.join(helper.fastq_path, '10015AT', '10015ATA0001L05', 'this.fastq.gz')
    expected = '--nogroup -t 4 -q ' + test_fastq

    assert commands.fastqc(test_fastq).endswith(expected)


def test_bcbio():
    run_yaml = os.path.join(helper.assets_path, 'run.yaml')
    expected = '%s -n 16 --workdir %s' % (run_yaml, helper.assets_path)

    assert commands.bcbio(run_yaml, helper.assets_path, cores=16).endswith(expected)
