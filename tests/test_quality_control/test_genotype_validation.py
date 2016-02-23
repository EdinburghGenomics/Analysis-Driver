from unittest.mock import patch, Mock, MagicMock
from tests.test_analysisdriver import TestAnalysisDriver
import os.path
import builtins
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

    @patch('analysis_driver.executor.execute')
    def test__vcf_validation(self, mocked_execute):
        vcf_file = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_expected_genotype.vcf')
        genotype_vcf = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        validation_results = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.txt')
        vc = self.validator.validation_cfg
        sample2genotype={self.sample_id:vcf_file}
        with patch('os.path.isfile', side_effect = [True, False]):
            self.validator._vcf_validation(sample2genotype)
        command_gatk = ' '.join(
                ['java -Xmx4G -jar %s' % vc.get('gatk'),
                 '-T GenotypeConcordance',
                 '-eval:VCF %s ' % vcf_file,
                 '-comp:VCF %s ' % genotype_vcf,
                 '-R %s' % vc.get('reference'),
                 ' > %s' % validation_results])
        command_index = '{tabix} -p vcf {vcf}'.format(tabix=vc.get('tabix'), vcf=genotype_vcf)
        #Call the index and the actuall validation
        assert mocked_execute.call_count == 2
        mocked_execute.assert_any_call([command_index], job_name='index_vcf', run_id=self.sample_id, cpus=1, mem=4)
        mocked_execute.assert_called_with([command_gatk], job_name='genotype_concordance', run_id=self.sample_id, cpus=4, mem=8, log_command=False)

        mocked_execute.reset_mock()
        with patch('os.path.isfile', side_effect = [True, True]):
            self.validator._vcf_validation(sample2genotype)
        #No call to the index generation and the actuall validation
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_with([command_gatk], job_name='genotype_concordance', run_id=self.sample_id, cpus=4, mem=8, log_command=False)

    def test__merge_validation_results(self):
        sample2genotype_validation={"T00001P001A01":os.path.join(self.assets_path, 'sample_data', "T00001P001A01-validation.txt"),
                                    "T00001P001A02":os.path.join(self.assets_path, 'sample_data', "T00001P001A02-validation.txt")}
        open_file1 = open(sample2genotype_validation.get("T00001P001A01"))
        open_file2 = open(sample2genotype_validation.get("T00001P001A02"))
        with patch.object(builtins, 'open') as mocked_open:
            mock_open_write = MagicMock()
            mocked_open.side_effect= [open_file1, open_file2, mock_open_write]
            self.validator._merge_validation_results(sample2genotype_validation)
            print(mocked_open.mock_calls)
            print(mock_open_write.mock_calls)

        open_file1.close()
        open_file2.close()

    @patch('analysis_driver.executor.execute')
    def test__rename_expected_genotype(self, mocked_execute):
        genotype_vcf = os.path.join('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        sample_name = 'test_sample'
        self.validator._rename_expected_genotype({sample_name:genotype_vcf})
        vc = self.validator.validation_cfg
        command = "{bcftools} reheader -s <(echo {sample_name}) {genotype_vcf} > {genotype_vcf}.tmp; mv {genotype_vcf}.tmp {genotype_vcf}"
        command = command.format(bcftools=vc.get('bcftools'), sample_name=sample_name, genotype_vcf=genotype_vcf)
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_once_with([command], job_name='genotype_rename', run_id=self.sample_id, cpus=4, mem=8, log_command=False)

    def test__genotype_validation(self):
        pass