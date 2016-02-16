__author__ = 'tcezard'
from unittest.mock import patch
from tests.test_analysisdriver import TestAnalysisDriver
import os.path
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import GenotypeValidation


class TestGenotypeValidation(TestAnalysisDriver):

    def setUp(self):
        self.sample_id = 'test_sample'
        self.fastq_files = [os.path.join('samples', self.sample_id, 'fastq_R1.fastq.gz'),
                       os.path.join('samples', self.sample_id, 'fastq_R2.fastq.gz')]
        self.validator = GenotypeValidation(fastq_files=self.fastq_files, sample_id=self.sample_id)

    def test_bwa_aln(self):
        cmd = self.validator._bwa_aln(
            self.fastq_files,
            self.sample_id,
            self.sample_id + '.bam',
            cfg['genotype-validation']['reference']
        )
        expected = (
            cfg['tools']['bwa'],
            'sampe -r \'@RG\\tID:1\\tSM:%s\'' % self.sample_id,
            cfg['genotype-validation']['reference'],
            '<(' + cfg['tools']['bwa'], 'aln', cfg['genotype-validation']['reference'], self.fastq_files[0] + ')',
            '<(' + cfg['tools']['bwa'], 'aln', cfg['genotype-validation']['reference'], self.fastq_files[1] + ')',
            self.fastq_files[0],
            self.fastq_files[1],
            '|',
            cfg['tools']['samblaster'], '--removeDups',
            '|',
            cfg['tools']['samtools'], 'view -F 4 -Sb -',
            '|',
            cfg['tools']['sambamba'], 'sort -t 16 -o ', self.sample_id + '.bam', '/dev/stdin'
        )
        assert cmd == ' '.join(expected)

    @patch('analysis_driver.quality_control.GenotypeValidation._bwa_aln', return_value='long_bwa_command')
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
        expected_vcf = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        command = ' '.join(
            ['java -Xmx4G -jar %s' % cfg.query('tools', 'gatk'),
             '-T UnifiedGenotyper',
             '-nt 4',
             '-R %s' % cfg.query('genotype-validation', 'reference'),
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

    @patch('analysis_driver.executor.execute')
    def test__vcf_validation(self, mocked_execute):
        vcf_file = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_expected_genotype.vcf')
        genotype_vcf = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        self.validator._vcf_validation(vcf_file, genotype_vcf)
        command = ' '.join(
            ['java -Xmx4G -jar %s' % cfg.query('tools', 'gatk'),
             '-T GenotypeConcordance',
             '-eval:VCF %s ' % vcf_file,
             '-comp:VCF %s ' % genotype_vcf,
             '-R %s' % cfg.query('genotype-validation', 'reference'),
             ' > %s' % self.validator.validation_results])
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_once_with([command], job_name='genotype_concordance', run_id=self.sample_id, cpus=4, mem=8, log_command=False)

    def test_genotype_validation(self):
        pass
