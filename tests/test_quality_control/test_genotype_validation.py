from os.path import join as pjoin
from unittest.mock import patch, mock_open, call
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import GenotypeValidation
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance


class TestGenotypeValidation(QCTester):
    def setUp(self):
        super().setUp()
        self.sample_id = 'test_sample'
        self.fastq_files = [
            pjoin('samples', self.sample_id, 'fastq_R1.fastq.gz'),
            pjoin('samples', self.sample_id, 'fastq_R2.fastq.gz')
        ]
        self.validator = GenotypeValidation(
            self.dataset,
            working_dir=pjoin(cfg['jobs_dir'], self.sample_id),
            fastq_files=self.fastq_files
        )

    def test_bwa_aln(self):
        reference = cfg['genotype-validation']['reference']
        cmd = self.validator._bwa_aln(
            self.fastq_files,
            self.sample_id,
            self.sample_id + '.bam',
            reference
        )
        expected = (
            cfg['tools']['bwa'],
            'sampe -r \'@RG\\tID:1\\tSM:%s\'' % self.sample_id,
            cfg['genotype-validation']['reference'],
            '<(' + cfg['tools']['bwa'], 'aln', reference, self.fastq_files[0] + ')',
            '<(' + cfg['tools']['bwa'], 'aln', reference, self.fastq_files[1] + ')',
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

    @patch('analysis_driver.dataset_scanner.rest_communication')
    @patch('analysis_driver.quality_control.GenotypeValidation._bwa_aln', return_value='long_bwa_command')
    @patch('analysis_driver.executor.execute')
    def test__bwa_alignment(self, mocked_execute, mocked_bwa_aln, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        expected_bam = self.validator._bwa_alignment()
        assert expected_bam == pjoin('path/to/jobs/', self.sample_id, self.sample_id + '_geno_val.bam')
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_once_with(
            'long_bwa_command',
            job_name='alignment_bwa',
            cpus=4,
            working_dir=pjoin(cfg['jobs_dir'], self.sample_id),
            mem=8
        )

    @patch('analysis_driver.dataset_scanner.rest_communication')
    @patch('analysis_driver.executor.execute')
    def test__snp_calling(self, mocked_execute, mocked_rest):
        bam_file = pjoin('path/to/bam', self.sample_id, self.sample_id + '_geno_val.bam')
        expected_vcf = pjoin('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
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
        mocked_execute.assert_called_once_with(
            command,
            job_name='snpcall_gatk',
            working_dir=pjoin(cfg['jobs_dir'], self.sample_id),
            cpus=4,
            mem=4
        )

    @patch('analysis_driver.dataset_scanner.rest_communication')
    @patch('analysis_driver.executor.execute')
    def test__vcf_validation(self, mocked_execute, mocked_rest):
        work_dir = pjoin(cfg.query('jobs_dir'), self.sample_id)
        vcf_file = pjoin(work_dir, self.sample_id + '_expected_genotype.vcf')
        genotype_vcf = pjoin(work_dir, self.sample_id + '_genotype_validation.vcf.gz')
        validation_results = pjoin(work_dir, self.sample_id + '_genotype_validation.txt')

        sample2genotype = {self.sample_id: vcf_file}
        with patch('os.path.isfile', side_effect=[True, False]):
            self.validator._vcf_validation(sample2genotype)
        command_gatk = ' '.join(
                ['java -Xmx4G -jar %s' % cfg['tools']['gatk'],
                 '-T GenotypeConcordance',
                 '-eval:VCF %s ' % genotype_vcf,
                 '-comp:VCF %s ' % vcf_file,
                 '-R %s' % cfg.query('genotype-validation', 'reference'),
                 ' > %s' % validation_results])
        command_index = '{tabix} -p vcf {vcf}'.format(tabix=cfg.query('tools', 'tabix'), vcf=genotype_vcf)
        # Call the index and the actuall validation
        assert mocked_execute.call_count == 2
        print(mocked_execute.call_args_list)
        mocked_execute.assert_any_call(command_index, job_name='index_vcf', working_dir=work_dir, cpus=1, mem=4)
        mocked_execute.assert_called_with(
            command_gatk,
            job_name='genotype_concordance',
            working_dir=work_dir,
            cpus=4,
            mem=8,
            log_commands=False
        )

        mocked_execute.reset_mock()
        with patch('os.path.isfile', side_effect=[True, True]):
            self.validator._vcf_validation(sample2genotype)
        # No call to the index generation and the actuall validation
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_with(
            command_gatk,
            job_name='genotype_concordance',
            working_dir=work_dir,
            cpus=4,
            mem=8,
            log_commands=False
        )

    def test__merge_validation_results(self):
        sample2genotype_validation = dict([
                ('T00001P001A01', pjoin(self.assets_path, 'sample_data', 'T00001P001A01-validation.txt')),
                ('T00001P001A02', pjoin(self.assets_path, 'sample_data', 'T00001P001A02-validation.txt'))
        ])
        o1 = parse_genotype_concordance(sample2genotype_validation.get('T00001P001A01'))
        o2 = parse_genotype_concordance(sample2genotype_validation.get('T00001P001A02'))

        with patch('analysis_driver.quality_control.genotype_validation.parse_genotype_concordance', side_effect=[o1, o2]):
            mock_open_write = mock_open()
            with patch('builtins.open', mock_open_write):
                self.validator._merge_validation_results(sample2genotype_validation)
                calls = []
                calls.append(call('path/to/jobs/test_sample/test_sample_genotype_validation.txt', 'w'))
                calls.append(call().__enter__())
                calls.append(call().write('#:GATKTable:GenotypeConcordance_Counts:Per-sample concordance tables: comparison counts\n'))
                calls.append(call().write('Sample\tNO_CALL_NO_CALL\tNO_CALL_HOM_REF\tNO_CALL_HET\tNO_CALL_HOM_VAR\t'
                                          'NO_CALL_UNAVAILABLE\tNO_CALL_MIXED\tHOM_REF_NO_CALL\tHOM_REF_HOM_REF\t'
                                          'HOM_REF_HET\tHOM_REF_HOM_VAR\tHOM_REF_UNAVAILABLE\tHOM_REF_MIXED\tHET_NO_CALL'
                                          '\tHET_HOM_REF\tHET_HET\tHET_HOM_VAR\tHET_UNAVAILABLE\tHET_MIXED\tHOM_VAR_NO_CALL'
                                          '\tHOM_VAR_HOM_REF\tHOM_VAR_HET\tHOM_VAR_HOM_VAR\tHOM_VAR_UNAVAILABLE\t'
                                          'HOM_VAR_MIXED\tUNAVAILABLE_NO_CALL\tUNAVAILABLE_HOM_REF\tUNAVAILABLE_HET\t'
                                          'UNAVAILABLE_HOM_VAR\tUNAVAILABLE_UNAVAILABLE\tUNAVAILABLE_MIXED\t'
                                          'MIXED_NO_CALL\tMIXED_HOM_REF\tMIXED_HET\tMIXED_HOM_VAR\tMIXED_UNAVAILABLE\t'
                                          'MIXED_MIXED\tMismatching_Alleles\n'))
                calls.append(call().write('T00001P001A01\t3\t0\t0\t0\t3859\t0\t0\t9\t0\t0\t34303\t0\t0\t0\t13\t0\t127\t'
                                          '0\t1\t0\t0\t6\t14\t0\t0\t0\t0\t0\t97\t0\t0\t0\t0\t0\t0\t0\t0\n'))
                calls.append(call().write('T00001P001A02\t3\t0\t0\t0\t3859\t0\t0\t9\t0\t0\t34303\t0\t0\t0\t13\t0\t127\t'
                                          '0\t1\t0\t0\t6\t14\t0\t0\t0\t0\t0\t97\t0\t0\t0\t0\t0\t0\t0\t0\n'))
                calls.append(call().__exit__(None, None, None))
                assert mock_open_write.mock_calls == calls

    @patch('analysis_driver.executor.execute')
    def test__rename_expected_genotype(self, mocked_execute):
        genotype_vcf = pjoin('path/to/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        sample_name = 'test_sample'
        work_dir = pjoin(cfg.query('jobs_dir'), sample_name)
        self.validator._rename_expected_genotype({sample_name:genotype_vcf})
        command = "{bcftools} reheader -s <(echo {sample_name}) {genotype_vcf} > {genotype_vcf}.tmp; mv {genotype_vcf}.tmp {genotype_vcf}"
        command = command.format(bcftools=cfg.query('tools','bcftools'), sample_name=sample_name, genotype_vcf=genotype_vcf)
        mocked_execute.assert_called_once_with(command, job_name='genotype_rename', working_dir=work_dir, cpus=4, mem=8, log_commands=False)

    def test_genotype_validation(self):
        pass
