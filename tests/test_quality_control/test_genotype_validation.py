from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver
import os.path
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import GenotypeValidation

__author__ = 'tcezard'


class TestGenotypeValidation(TestAnalysisDriver):

    def setUp(self):
        self.sample_id = 'test_sample'
        self.fastq_files = [os.path.join('samples', self.sample_id, 'fastq_R1.fastq.gz'),
                       os.path.join('samples', self.sample_id, 'fastq_R2.fastq.gz')]
        self.validator = GenotypeValidation(fastqs_files=self.fastq_files, sample_id=self.sample_id)

    def test_bwa_aln(self):
        vc = self.validator.validation_cfg
        cmd = self.validator._bwa_aln(
                self.fastq_files,
                self.sample_id,
                self.sample_id + '.bam',
                vc['reference']
            )
        expected = (
            vc['bwa'],
            'sampe -r \'@RG\\tID:1\\tSM:%s\'' % self.sample_id,
            vc['reference'],
            '<(' + vc['bwa'], 'aln', vc['reference'], self.fastq_files[0] + ')',
            '<(' + vc['bwa'], 'aln', vc['reference'], self.fastq_files[1] + ')',
            self.fastq_files[0],
            self.fastq_files[1],
            '|',
            vc['samblaster'], '--removeDups',
            '|',
            vc['samtools'], 'view -F 4 -Sb -',
            '|',
            vc['sambamba'], 'sort -t 16 -o ', self.sample_id + '.bam', '/dev/stdin'
        )
        assert cmd == ' '.join(expected)

    @patch('analysis_driver.quality_control.GenotypeValidation._bwa_aln', return_value = 'long_bwa_command')
    @patch('analysis_driver.executor.execute')
    def test__bwa_alignment(self, mocked_execute, mocked_bwa_aln):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        expected_bam = self.validator._bwa_alignment()
        assert expected_bam == os.path.join('path/to/jobs/', self.sample_id, self.sample_id + '_geno_val.bam')
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_once_with(['long_bwa_command'], job_name='alignment_bwa', cpus=4,
                                               run_id='test_sample', mem=8)

    @patch('analysis_driver.executor.execute')
    def test__snp_calling(self, mocked_execute):
        bam_file = os.path.join('path/to/bam', self.sample_id, self.sample_id + '_geno_val.bam')
        vc = self.validator.validation_cfg
        expected_vcf = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        command = ' '.join(
                ['java -Xmx4G -jar %s' % vc.get('gatk'),
                 '-T UnifiedGenotyper',
                 '-nt 4',
                 '-R %s' % vc.get('reference'),
                 ' --standard_min_confidence_threshold_for_calling 30.0',
                 '--standard_min_confidence_threshold_for_emitting 0',
                 '-out_mode EMIT_ALL_SITES',
                 '-I %s' % bam_file,
                 '-o %s' % expected_vcf]
        )
        output_vcf = self.validator._snp_calling(bam_file)
        assert output_vcf == expected_vcf
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_once_with([command], job_name='snpcall_gatk', run_id=self.sample_id, cpus=4, mem=4)

    def test__vcf_validation(self):
        pass

    def test__genotype_validation(self):
        pass