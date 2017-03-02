import shutil
import os
from unittest.mock import Mock, patch
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import transfer_data, util
from analysis_driver.config import default as cfg


def patched_get_user_sample_name(sample_id):
    return patch('egcg_core.clarity.get_user_sample_name', return_value=sample_id)


def patched_find_project_from_sample(sample_id):
    return patch('egcg_core.clarity.find_project_name_from_sample', return_value='proj_' + sample_id)


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
        self._create_pseudo_links()

    def tearDown(self):
        for p in self.param_remappings:
            cfg.content[p['name']] = p['original']

        shutil.rmtree(self._to_dir)
        shutil.rmtree(self._pseudo_links)

    @property
    def _to_dir(self):
        return os.path.join(self.data_output, 'to', '')

    def _create_pseudo_links(self):
        os.makedirs(self._pseudo_links, exist_ok=True)
        for f in [
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
        ]:
            open(os.path.join(self._pseudo_links, f), 'a').close()

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


    @property
    def _pseudo_links(self):
        return os.path.join(self.data_output, 'pseudo_links')

    def test_output_sample_data(self):
        with patched_find_project_from_sample(self.sample_id), \
                patch('analysis_driver.transfer_data.archive_management.archive_directory', return_value=True):
            exit_status = transfer_data.output_sample_data(
                sample_id=self.sample_id,
                source_dir=self._pseudo_links,
                output_dir=self._to_dir
            )
        output_files = os.path.join(self._to_dir, 'proj_' + self.sample_id, self.sample_id)

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
        expected = (
            'bcbio_prepare_samples.py --out ' +
            os.path.join(self.assets_path, 'merged') +
            ' --csv ' +
            bcbio_csv
        )
        assert expected in cmd


class TestBCLValidator(TestAnalysisDriver):
    def setUp(self):
        run_info = Mock(
            tiles=('s_1_1101', 's_2_1101', 's_1_1102', 's_2_1102'),
            mask=Mock(reads=[Mock(attrib={'NumCycles': '3'})])
        )
        self.job_dir = os.path.join(TestAnalysisDriver.assets_path, 'bcl_validation')
        validation_log = os.path.join(self.job_dir, 'bcl_validation.log')
        if os.path.isfile(validation_log):
            os.remove(validation_log)

        self.val = util.BCLValidator(self.job_dir, run_info, validation_log)

    @patch('analysis_driver.util.BCLValidator._all_cycles_from_interop')
    def test_get_bcl_files_to_check(self, mocked_cycles):
        mocked_cycles.return_value = [1, 1, 1]  # no completed cycles
        assert self.val.get_bcls_to_check() == []

        mocked_cycles.return_value.extend([1, 2, 2, 2])  # completed cycle 1, but not cycle 2
        assert self.val.get_bcls_to_check() == [
            os.path.join(self.val.basecalls_dir, f)
            for f in ('L001/C1.1/s_1_1101.bcl.gz', 'L002/C1.1/s_2_1101.bcl.gz',
                      'L001/C1.1/s_1_1102.bcl.gz', 'L002/C1.1/s_2_1102.bcl.gz')
        ]

        mocked_cycles.return_value.extend([2, 3, 3, 3, 3])
        obs = self.val.get_bcls_to_check()
        exp = [
            os.path.join(self.val.basecalls_dir, f)
            for f in (
                'L001/C1.1/s_1_1101.bcl.gz', 'L001/C1.1/s_1_1102.bcl.gz',
                'L001/C2.1/s_1_1101.bcl.gz', 'L001/C2.1/s_1_1102.bcl.gz',
                'L001/C3.1/s_1_1101.bcl.gz', 'L001/C3.1/s_1_1102.bcl.gz',
                'L002/C1.1/s_2_1101.bcl.gz', 'L002/C1.1/s_2_1102.bcl.gz',
                'L002/C2.1/s_2_1101.bcl.gz', 'L002/C2.1/s_2_1102.bcl.gz',
                'L002/C3.1/s_2_1101.bcl.gz', 'L002/C3.1/s_2_1102.bcl.gz'
            )
        ]
        assert sorted(obs) == sorted(exp)

    @patch('analysis_driver.util.executor.execute', return_value=Mock(join=Mock(return_value=0)))
    def test_run_bcl_check(self, mocked_execute):
        with patch('analysis_driver.util.BCLValidator._all_cycles_from_interop',
                   return_value=[1, 1, 1, 1]):
            bcls = self.val.get_bcls_to_check()
        self.val.run_bcl_check(bcls, self.job_dir, slice_size=2)

        mocked_execute.assert_called_with(
            '\n'.join('check_bcl ' + os.path.join(self.job_dir, f) for f in bcls[:2]),
            '\n'.join('check_bcl ' + os.path.join(self.job_dir, f) for f in bcls[2:]),
            prelim_cmds=[self.val.validate_expr],
            job_name='bcl_validation',
            working_dir=self.job_dir,
            log_commands=False,
            cpus=1,
            mem=6
        )

        e = util.executor.SlurmExecutor(
            '\n'.join('check_bcl ' + os.path.join(self.job_dir, f) for f in bcls[:2]),
            '\n'.join('check_bcl ' + os.path.join(self.job_dir, f) for f in bcls[2:]),
            prelim_cmds=[self.val.validate_expr],
            job_name='bcl_validation',
            working_dir=self.job_dir,
            log_commands=False,
            cpus=1,
            mem=6
        )
        e.write_script()

    def test_run_bcl_check_local(self):
        with patch('analysis_driver.util.BCLValidator._all_cycles_from_interop',
                   return_value=[1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]):
            bcls = self.val.get_bcls_to_check()
        self.val.run_bcl_check_local(bcls)
        assert self.val.read_invalid_files() == [
            os.path.join(self.val.basecalls_dir, 'L002', 'C3.1', 's_2_1101.bcl.gz'),
            os.path.join(self.val.basecalls_dir, 'L001', 'C3.1', 's_1_1102.bcl.gz')
        ]

    def test_cycles_from_interop(self):
        interop_dir = os.path.join(self.job_dir, 'InterOp')
        os.makedirs(interop_dir, exist_ok=True)
        assert self.val._all_cycles_from_interop(self.job_dir) == []  # no ExtractionMetrics
        open(os.path.join(interop_dir, 'ExtractionMetricsOut.bin'), 'w').close()
        assert self.val._all_cycles_from_interop(self.job_dir) == []  # empty ExtractionMetrics
