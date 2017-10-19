from os.path import join
from unittest.mock import patch, mock_open, call
from analysis_driver.util.bash_commands import export_env_vars
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.config import default as cfg
from analysis_driver.quality_control.genotype_validation import GenotypeValidation
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance

patched_execute = patch('egcg_core.executor.execute')
ppath = 'analysis_driver.quality_control.genotype_validation.'


class TestGenotypeValidation(QCTester):
    def setUp(self):
        super().setUp()
        self.fastq_files = [
            join('samples', self.sample_id, 'fastq_R1.fastq.gz'),
            join('samples', self.sample_id, 'fastq_R2.fastq.gz')
        ]
        self.validator = GenotypeValidation(
            dataset=self.dataset,
            fq_pattern=join('samples', self.sample_id, 'fastq_R?.fastq.gz')
        )

    def test_bwa_aln(self):
        observed = self.validator._bwa_aln(self.fastq_files)
        expected = (
            "path/to/bwa sampe -r '@RG\\tID:1\\tSM:test_sample' path/to/32_snps_600bp.fa "
            '<(path/to/bwa aln path/to/32_snps_600bp.fa {fq1}) '
            '<(path/to/bwa aln path/to/32_snps_600bp.fa {fq2}) '
            '{fq1} {fq2} | path/to/samblaster --removeDups | path/to/samtools view -F 4 -Sb - | '
            'path/to/sambamba sort -t 16 -o tests/assets/jobs/test_sample/test_sample_geno_val.bam /dev/stdin'
        ).format(fq1=self.fastq_files[0], fq2=self.fastq_files[1])
        assert observed == expected

    @patch(ppath + 'GenotypeValidation._bwa_aln', return_value='long_bwa_command')
    @patched_execute
    def test_bwa_alignment(self, mocked_execute, mocked_bwa_aln):
        self.validator._bwa_alignment()
        assert self.validator.output_bam == join('tests/assets/jobs/test_sample/test_sample_geno_val.bam')
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_once_with(
            'long_bwa_command',
            job_name='alignment_bwa',
            cpus=4,
            working_dir=join(cfg['jobs_dir'], self.sample_id),
            mem=8
        )

    @patched_execute
    def test_snp_calling(self, mocked_execute):
        bam_file = join('path/to/bam', self.sample_id, self.sample_id + '_geno_val.bam')
        expected_vcf = join('tests/assets/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        self.validator._snp_calling(bam_file)
        assert self.validator.seq_vcf_file == expected_vcf
        command = (
            'java -Djava.io.tmpdir={work_dir} -XX:+UseSerialGC -Xmx4G -jar path/to/gatk '
            '-T UnifiedGenotyper -nt 4 -R path/to/32_snps_600bp.fa '
            '--standard_min_confidence_threshold_for_calling 30.0 '
            '--standard_min_confidence_threshold_for_emitting 0 -out_mode EMIT_ALL_SITES '
            '-I {bam_file} -o {vcf_file}'
        ).format(work_dir=join(cfg['jobs_dir'], self.sample_id), bam_file=bam_file, vcf_file=expected_vcf)

        mocked_execute.assert_called_with(
            command,
            prelim_cmds=export_env_vars(),
            job_name='snpcall_gatk',
            working_dir=join(cfg['jobs_dir'], self.sample_id),
            cpus=4,
            mem=4
        )

    @patched_execute
    def test_vcf_validation(self, mocked_execute):
        work_dir = self.validator.job_dir
        vcf_file = join(work_dir, 'test_sample_expected_genotype.vcf')
        genotype_vcf = join(work_dir, 'test_sample_genotype_validation.vcf.gz')
        validation_results = join(work_dir, 'test_sample_genotype_validation.txt')

        sample2genotype = {self.sample_id: vcf_file}
        with patch('os.path.isfile', side_effect=[True, False]):
            assert self.validator._vcf_validation(sample2genotype) == {
                'test_sample': 'tests/assets/jobs/test_sample/test_sample_genotype_validation.txt'
            }

        command_gatk = (
            'java -Djava.io.tmpdir={work_dir} -XX:+UseSerialGC -Xmx4G -jar path/to/gatk '
            '-T GenotypeConcordance -eval:VCF {genovcf} -comp:VCF {sample_vcf} '
            '-R path/to/32_snps_600bp.fa > {val_results}'
        ).format(work_dir=work_dir, genovcf=genotype_vcf, sample_vcf=vcf_file, val_results=validation_results)
        # Call the index and the actual validation
        assert mocked_execute.call_count == 2
        mocked_execute.assert_any_call(
            'path/to/tabix -p vcf ' + genotype_vcf,
            job_name='index_vcf',
            working_dir=work_dir,
            cpus=1,
            mem=4
        )
        mocked_execute.assert_called_with(
            command_gatk,
            prelim_cmds=export_env_vars(),
            job_name='genotype_concordance',
            working_dir=work_dir,
            cpus=4,
            mem=8,
            log_commands=False
        )

        mocked_execute.reset_mock()
        with patch('os.path.isfile', side_effect=[True, True]):
            self.validator._vcf_validation(sample2genotype)
        # No call to the index generation and the actual validation
        assert mocked_execute.call_count == 1
        mocked_execute.assert_called_with(
            command_gatk,
            prelim_cmds=export_env_vars(),
            job_name='genotype_concordance',
            working_dir=work_dir,
            cpus=4,
            mem=8,
            log_commands=False
        )

    def test_merge_validation_results(self):
        sample2genotype_validation = {
            'T00001P001A01': join(self.assets_path, 'sample_data', 'T00001P001A01-validation.txt'),
            'T00001P001A02': join(self.assets_path, 'sample_data', 'T00001P001A02-validation.txt')
        }
        o1 = parse_genotype_concordance(sample2genotype_validation['T00001P001A01'])
        o2 = parse_genotype_concordance(sample2genotype_validation['T00001P001A02'])

        with patch(ppath + 'parse_genotype_concordance', side_effect=[o1, o2]):
            mock_open_write = mock_open()
            with patch('builtins.open', mock_open_write):
                self.validator._merge_validation_results(sample2genotype_validation)
                calls = [
                    call('tests/assets/jobs/test_sample/test_sample_genotype_validation.txt', 'w'),
                    call().__enter__(),
                    call().write('#:GATKTable:GenotypeConcordance_Counts:Per-sample concordance tables: comparison counts\n'),
                    call().write(
                        'Sample\tNO_CALL_NO_CALL\tNO_CALL_HOM_REF\tNO_CALL_HET\tNO_CALL_HOM_VAR\t'
                        'NO_CALL_UNAVAILABLE\tNO_CALL_MIXED\tHOM_REF_NO_CALL\tHOM_REF_HOM_REF\t'
                        'HOM_REF_HET\tHOM_REF_HOM_VAR\tHOM_REF_UNAVAILABLE\tHOM_REF_MIXED\tHET_NO_CALL'
                        '\tHET_HOM_REF\tHET_HET\tHET_HOM_VAR\tHET_UNAVAILABLE\tHET_MIXED\tHOM_VAR_NO_CALL'
                        '\tHOM_VAR_HOM_REF\tHOM_VAR_HET\tHOM_VAR_HOM_VAR\tHOM_VAR_UNAVAILABLE\t'
                        'HOM_VAR_MIXED\tUNAVAILABLE_NO_CALL\tUNAVAILABLE_HOM_REF\tUNAVAILABLE_HET\t'
                        'UNAVAILABLE_HOM_VAR\tUNAVAILABLE_UNAVAILABLE\tUNAVAILABLE_MIXED\t'
                        'MIXED_NO_CALL\tMIXED_HOM_REF\tMIXED_HET\tMIXED_HOM_VAR\tMIXED_UNAVAILABLE\t'
                        'MIXED_MIXED\tMismatching_Alleles\n'
                    ),
                    call().write(
                        'T00001P001A01\t3\t0\t0\t0\t3859\t0\t0\t9\t0\t0\t34303\t0\t0\t0\t13\t0\t127\t'
                        '0\t1\t0\t0\t6\t14\t0\t0\t0\t0\t0\t97\t0\t0\t0\t0\t0\t0\t0\t0\n'
                    ),
                    call().write(
                        'T00001P001A02\t3\t0\t0\t0\t3859\t0\t0\t9\t0\t0\t34303\t0\t0\t0\t13\t0\t127\t'
                        '0\t1\t0\t0\t6\t14\t0\t0\t0\t0\t0\t97\t0\t0\t0\t0\t0\t0\t0\t0\n'
                    ),
                    call().__exit__(None, None, None),
                ]
                assert mock_open_write.mock_calls == calls

    @patched_execute
    def test_rename_expected_genotype(self, mocked_execute):
        genotype_vcf = join('tests/assets/jobs', self.sample_id, self.sample_id + '_genotype_validation.vcf.gz')
        self.validator._rename_expected_genotype({self.sample_id: genotype_vcf})
        command = ('path/to/bcftools reheader -s <(echo test_sample) {genotype_vcf} > {genotype_vcf}.tmp; '
                   'mv {genotype_vcf}.tmp {genotype_vcf}').format(genotype_vcf=genotype_vcf)
        mocked_execute.assert_called_once_with(
            command,
            job_name='genotype_rename',
            working_dir=self.validator.job_dir,
            cpus=4,
            mem=8,
            log_commands=False
        )

    def test_run(self):
        self.validator.check_neighbour = True
        self.validator.check_project = True
        self.validator.list_samples = ['another_sample']

        with patch(ppath + 'GenotypeValidation._bwa_alignment') as bwa,\
            patch(ppath + 'GenotypeValidation._snp_calling') as snp,\
            patch(ppath + 'clarity.get_samples_arrived_with', return_value=['a_sample']),\
            patch(ppath + 'clarity.get_samples_genotyped_with', return_value=['a_sample']),\
            patch(ppath + 'clarity.get_samples_sequenced_with', return_value=['a_sample']),\
            patch(ppath + 'clarity.find_project_name_from_sample'),\
            patch(ppath + 'clarity.get_sample_names_from_project', return_value=['a_sample']),\
            patch(ppath + 'clarity.get_sample_genotype', return_value='test.vcf'),\
            patch(ppath + 'GenotypeValidation._rename_expected_genotype') as rename,\
            patch(ppath + 'GenotypeValidation._vcf_validation', return_value=[None, None]) as vcf_val,\
            patch(ppath + 'GenotypeValidation._merge_validation_results') as merge:
                self.validator._run()

        assert bwa.call_count == 1
        snp.assert_called_with('tests/assets/jobs/test_sample/test_sample_geno_val.bam')

        fake_sample2genotype = {k: 'test.vcf' for k in ('test_sample', 'a_sample', 'another_sample')}
        rename.assert_called_with(fake_sample2genotype)
        vcf_val.assert_called_with(fake_sample2genotype)

        merge.assert_called_with([None, None])
        assert self.validator.validation_results == [None, None]
