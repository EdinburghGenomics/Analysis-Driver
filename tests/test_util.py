from analysis_driver.util.bash_commands import sickle_paired_end_in_place

__author__ = 'mwham'
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import util, transfer_data
from analysis_driver.config import default as cfg
import shutil
import os.path


def patched_get_user_sample_name(sample_id):
    return patch('analysis_driver.clarity.get_user_sample_name', return_value=sample_id)


def patched_find_project_from_sample(sample_id):
    return patch('analysis_driver.clarity.find_project_from_sample', return_value='proj_' + sample_id)


class TestUtil(TestAnalysisDriver):
    def test_find_fastqs(self):
        fastqs = util.find_fastqs(self.fastq_path, '10015AT', '10015AT0001')
        for file_name in ['10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz',
                          '10015AT0001_S6_L005_R1_001.fastq.gz', '10015AT0001_S6_L005_R2_001.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs

    def test_find_fastqs_with_lane(self):
        fastqs = util.find_fastqs(self.fastq_path, '10015AT', '10015AT0001', lane=4)
        for file_name in ['10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs

    def test_find_all_fastqs(self):
        fastqs = util.find_all_fastqs(self.fastq_path)
        for file_name in ['10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz']:
            assert os.path.join(
                self.fastq_path, '10015AT', '10015AT0001', file_name
            ) in fastqs


    def find_all_fastq_pairs(self):
        fastqs = util.find_all_fastq_pairs(self.fastq_path)
        for f1, f2 in [('10015AT0001_S6_L004_R1_001.fastq.gz', '10015AT0001_S6_L004_R2_001.fastq.gz'),
                       ('10015AT0001_S6_L005_R1_001.fastq.gz', '10015AT0001_S6_L005_R2_001.fastq.gz')]:
            fp1 = os.path.join(self.fastq_path, '10015AT', '10015AT0001', f1 )
            fp2 = os.path.join(self.fastq_path, '10015AT', '10015AT0001', f2 )
            assert (fp1, fp2) in fastqs

    def test_sickle_paired_end_in_place(self):
        expected_command = "path/to/sickle pe -f fastqfile1_R1.fastq.gz -r fastqfile1_R2.fastq.gz " \
                           "-o fastqfile1_R1.fastq_sickle.gz -p fastqfile1_R2.fastq_sickle.gz -s " \
                           "fastqfile1_R1.fastq_sickle_single.gz -q 5  -l 36  -x  -g -t sanger\n"\
                           "EXIT_CODE=$?\n$EXIT_CODE && mv fastqfile1_R1.fastq_sickle.gz fastqfile1_R1.fastq.gz\n" \
                           "$EXIT_CODE && mv fastqfile1_R2.fastq_sickle.gz fastqfile1_R2.fastq.gz\n" \
                           "$EXIT_CODE && rm fastqfile1_R1.fastq_sickle_single.gz\n(exit $EXIT_CODE)"
        cmd = sickle_paired_end_in_place(('fastqfile1_R1.fastq.gz', 'fastqfile1_R2.fastq.gz'))
        assert cmd == expected_command

class TestTransferData(TestAnalysisDriver):
    sample_id = '10015AT0001'

    def setUp(self):
        self.param_remappings = (
            {'name': 'output_dir', 'new': self._to_dir},
            {'name': 'jobs_dir', 'new': os.path.join(self.data_output, 'jobs')}
        )
        for p in self.param_remappings:
            p['original'] = cfg.get(p['name'])
            cfg.content[p['name']] = p['new']

        os.makedirs(self._to_dir, exist_ok=True)

    def tearDown(self):
        for p in self.param_remappings:
            cfg.content[p['name']] = p['original']

        shutil.rmtree(self._to_dir)

    @property
    def _to_dir(self):
        return os.path.join(self.data_output, 'to', '')

    def test_create_links(self):
        if not os.path.exists(os.path.join(self.data_output, 'jobs', self.sample_id)):
            os.makedirs(os.path.join(self.data_output, 'jobs', self.sample_id))

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

        dir_with_linked_files = os.path.join(self.data_output, 'linked_output_files')
        if os.path.isdir(dir_with_linked_files):
            shutil.rmtree(dir_with_linked_files)
        os.makedirs(dir_with_linked_files)

        with patched_get_user_sample_name(self.sample_id):
            list_of_linked_files = transfer_data.create_links_from_bcbio(
                self.sample_id,
                self.data_output,
                records,
                dir_with_linked_files
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

    def patched_rsync(self):
        return patch(
            'analysis_driver.transfer_data.rsync_from_to',
            return_value='rsync -rLD --size-only %s/ %s' % (  # the trailing slash is important...
                self._pseudo_links,
                os.path.join(self._to_dir, 'proj_' + self.sample_id, self.sample_id)
            )
        )

    @property
    def _pseudo_links(self):
        return os.path.join(self.data_output, 'pseudo_links')

    def test_output_sample_data(self):
        with patched_find_project_from_sample(self.sample_id), self.patched_rsync():
            exit_status = transfer_data.output_sample_data(
                sample_id=self.sample_id,
                source_dir=self._pseudo_links,
                output_dir=self._to_dir
            )
        output_files = os.path.join(
            self._to_dir,
            'proj_' + self.sample_id,
            self.sample_id
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

        o = sorted(os.listdir(output_files))
        assert exit_status == 0
        assert o == expected_outputs

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
