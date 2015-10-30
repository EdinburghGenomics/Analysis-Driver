__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import util
import shutil
from analysis_driver.driver import _output_data
from analysis_driver.reader import SampleSheet
import os.path

helper = TestAnalysisDriver()


def test_setup_bcbio_run():
    print('Setup_bcbio_run currently untestable')
    assert True


class TestFastqHandler(TestAnalysisDriver):
    def test_find_fastqs(self):
        fastqs = util.fastq_handler.find_fastqs(self.fastq_path, '10015AT', '10015AT0001')
        for file_name in ['10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz',
                          '10015AT0001_S6_L005_R1_001.fastq.gz', '10015AT0001_S6_L005_R2_001.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs

    def test_find_fastqs_with_lane(self):
        fastqs = util.fastq_handler.find_fastqs(self.fastq_path, '10015AT', '10015AT0001', lane=4)
        for file_name in ['10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs

    def test_find_all_fastqs(self):
        fastqs = util.fastq_handler.find_all_fastqs(self.fastq_path)
        for file_name in ['10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs


def test_output_data():
    sample_project = '10015AT'
    sample_id = '10015AT0001'
    destination = os.path.join(helper.data_output, 'output_data', sample_project)
    if not os.path.isdir(destination):
        os.makedirs(destination)
    for f in os.listdir(destination):
        shutil.rmtree(os.path.join(destination, f))
    assert not os.listdir(destination)

    records = [
        {
            'location': ['samples_{runfolder}-merged', 'final', '{sample_id}'],
            'basename': '{sample_id}-gatk-haplotype.vcf.gz',
            'new_name': '{sample_id}.g.vcf.gz'
        },
        {
            'location': ['samples_{runfolder}-merged', 'final', '{sample_id}'],
            'basename': '{sample_id}-gatk-haplotype.vcf.gz.tbi',
            'new_name': '{sample_id}.g.vcf.gz.tbi'
        },
        {
            'location': ['samples_{runfolder}-merged', 'final', '{sample_id}'],
            'basename': '{sample_id}-ready.bam',
            'new_name': '{sample_id}.bam'
        },
        {
            'location': ['samples_{runfolder}-merged', 'final', '{sample_id}'],
            'basename': '{sample_id}-ready.bam.bai',

            'new_name': '{sample_id}.bam.bai'
        },
        {
            'location': ['samples_{runfolder}-merged', 'final', '{sample_id}', 'qc', 'bamtools'],
            'basename': 'bamtools_stats.txt'
        },
        {
            'location': ['samples_{runfolder}-merged', 'work', 'align', '{sample_id}'],
            'basename': '*{sample_id}*-sort-highdepth-stats.yaml'
        },
        {
            'location': ['samples_{runfolder}-merged', 'work', 'align', '{sample_id}'],
            'basename': '*{sample_id}*-sort-callable.bed'
        },

        {'location': ['merged'], 'basename': '{sample_id}_R1.fastq.gz'},
        {'location': ['merged'], 'basename': '{sample_id}_R2.fastq.gz'},
        {'location': ['fastq', 'Stats'], 'basename': 'ConversionStats.xml'}
    ]

    sample_sheet = SampleSheet(helper.assets_path)
    exit_status = _output_data(
        sample_sheet,
        helper.data_output,
        os.path.join(helper.data_output, 'output_data'),
        records,
        query_lims=False
    )

    output_files = os.listdir(os.path.join(destination, sample_id))
    output_files.sort()

    expected_outputs = [
        '10015AT0001.bam',
        '10015AT0001.bam.bai',
        '10015AT0001.bam.bai.md5',
        '10015AT0001.bam.md5',
        '10015AT0001.g.vcf.gz',
        '10015AT0001.g.vcf.gz.md5',
        '10015AT0001.g.vcf.gz.tbi',
        '10015AT0001.g.vcf.gz.tbi.md5',
        '10015AT0001_R1.fastq.gz',
        '10015AT0001_R1.fastq.gz.md5',
        '10015AT0001_R2.fastq.gz',
        '10015AT0001_R2.fastq.gz.md5',
        '1_2015-10-16_samples_10015AT0001-merged-sort-callable.bed',
        '1_2015-10-16_samples_10015AT0001-merged-sort-callable.bed.md5',
        '1_2015-10-16_samples_10015AT0001-merged-sort-highdepth-stats.yaml',
        '1_2015-10-16_samples_10015AT0001-merged-sort-highdepth-stats.yaml.md5',
        'ConversionStats.xml',
        'ConversionStats.xml.md5',
        'bamtools_stats.txt',
        'bamtools_stats.txt.md5',
        'run_config.yaml'
    ]

    print(output_files)
    print(expected_outputs)
    assert exit_status == 0
    assert output_files == expected_outputs


def _join(*parts):
        return ''.join(parts)


def test_prep_samples_cmd():
    cmd = util.bcbio_prepare_samples_cmd(
        helper.assets_path,
        'a_sample_id',
        ['test_R1.fastq', 'test_R2.fastq'],
        user_sample_id='a_user_sample_id'
    )

    bcbio_csv = os.path.join(helper.assets_path, 'samples_a_sample_id.csv')
    with open(bcbio_csv) as f:
        content = f.read()
        print(content)
        assert content == (
            'samplename,description\n'
            'test_R1.fastq,a_user_sample_id\n'
            'test_R2.fastq,a_user_sample_id\n'
        )
    os.remove(bcbio_csv)
    assert not os.path.isfile(bcbio_csv)
    expected = _join(
        'bcbio_prepare_samples.py --out ',
        os.path.join(helper.assets_path, 'merged'),
        ' --csv ',
        bcbio_csv
    )
    assert expected in ' '.join(cmd)
