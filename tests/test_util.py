__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import util, transfer_data
from analysis_driver.config import default as cfg
import shutil
import os.path


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


class TestOutputData(TestAnalysisDriver):
    def setUp(self):
        self.param_remappings = (
            {'name': 'output_dir', 'new': os.path.join(self.data_output, 'to')},
            {'name': 'jobs_dir', 'new': os.path.join(self.data_output, 'jobs')}
        )

        for p in self.param_remappings:
            p['original'] = cfg.get(p['name'])
            cfg.content[p['name']] = p['new']

    def tearDown(self):
        for p in self.param_remappings:
            cfg.content[p['name']] = p['original']

    def test_create_links(self):
        sample_id = '10015AT0001'
        destination = os.path.join(self.data_output, 'to')
        if not os.path.exists(os.path.join(self.data_output, 'jobs', sample_id)):
            os.makedirs(os.path.join(self.data_output, 'jobs', sample_id))

        if not os.path.isdir(destination):
            os.makedirs(destination)

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

        if os.path.isdir(os.path.join(self.data_output, 'linked_output_files')):
            shutil.rmtree(os.path.join(self.data_output, 'linked_output_files'))

        list_of_linked_files = transfer_data.create_links_from_bcbio(
            sample_id,
            self.data_output,
            records,
            os.path.join(self.data_output, 'linked_output_files'),
            query_lims=False,
        )

        output_files = os.path.join(self.data_output, 'linked_output_files')

        expected_outputs = [
            '10015AT0001.bam',
            '10015AT0001.bam.bai',
            '10015AT0001.g.vcf.gz',
            '10015AT0001.g.vcf.gz.tbi',
            '10015AT0001_R1.fastq.gz',
            '10015AT0001_R2.fastq.gz',
            '1_2015-10-16_samples_10015AT0001-merged-sort-callable.bed',
            '1_2015-10-16_samples_10015AT0001-merged-sort-highdepth-stats.yaml',
            'ConversionStats.xml',
            'bamtools_stats.txt'
            # 'run_config.yaml'
        ]
        o = list(sorted(os.listdir(output_files)))
        assert len(list_of_linked_files) == len(expected_outputs)
        assert o == expected_outputs
        shutil.rmtree(output_files)
        assert not os.path.exists(output_files)

    def test_output_sample_data(self):
        sample_id = '10015AT0001'
        destination = os.path.join(self.data_output, 'to')
        if not os.path.isdir(destination):
            os.makedirs(destination)
        source = os.path.join(self.data_output, 'pseudo_links')

        exit_status = transfer_data.output_sample_data(
            sample_id,
            source,
            destination,
            query_lims=False,
            rsync_append=False
        )
        output_files = os.path.join(
            destination,
            'proj_' + sample_id,
            sample_id,
        )

        expected_outputs = [
            '10015AT0001.bam',
            '10015AT0001.bam.bai',
            '10015AT0001.g.vcf.gz',
            '10015AT0001.g.vcf.gz.tbi',
            '10015AT0001_R1.fastq.gz',
            '10015AT0001_R2.fastq.gz',
            '1_2015-10-16_samples_10015AT0001-merged-sort-callable.bed',
            '1_2015-10-16_samples_10015AT0001-merged-sort-highdepth-stats.yaml',
            'ConversionStats.xml',
            'bamtools_stats.txt'
            # 'run_config.yaml'
        ]

        o = list(sorted(os.listdir(output_files)))
        assert exit_status == 0
        assert o == expected_outputs
        shutil.rmtree(output_files)
        assert not os.path.exists(output_files)

    @staticmethod
    def _join(*parts):
            return ''.join(parts)

    def test_prep_samples_cmd(self):
        cmd = util.bcbio_prepare_samples_cmd(
            self.assets_path,
            'a_sample_id',
            ['test_R1.fastq', 'test_R2.fastq'],
            user_sample_id='a_user_sample_id'
        )

        bcbio_csv = os.path.join(self.assets_path, 'samples_a_sample_id.csv')
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
        expected = self._join(
            'bcbio_prepare_samples.py --out ',
            os.path.join(self.assets_path, 'merged'),
            ' --csv ',
            bcbio_csv
        )
        assert expected in cmd
