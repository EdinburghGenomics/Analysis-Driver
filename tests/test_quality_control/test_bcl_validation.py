import os
from unittest.mock import Mock, patch
from tests.test_analysisdriver import TestAnalysisDriver
from egcg_core import executor
from analysis_driver.quality_control import BCLValidator


class TestBCLValidator(TestAnalysisDriver):
    def setUp(self):
        run_info = Mock(
            tiles=('s_1_1101', 's_2_1101', 's_1_1102', 's_2_1102'),
            reads=[Mock(attrib={'NumCycles': '3'})]
        )
        self.job_dir = os.path.join(TestAnalysisDriver.assets_path, 'bcl_validation')
        validation_log = os.path.join(self.job_dir, 'bcl_validation.log')
        if os.path.isfile(validation_log):
            os.remove(validation_log)

        self.val = BCLValidator(self.job_dir, run_info, validation_log, Mock())

    @patch('analysis_driver.quality_control.BCLValidator._all_cycles_from_interop')
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

    @patch('analysis_driver.quality_control.bcl_validation.executor.execute', return_value=Mock(join=Mock(return_value=0)))
    def test_run_bcl_check(self, mocked_execute):
        with patch('analysis_driver.quality_control.BCLValidator._all_cycles_from_interop',
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

        e = executor.SlurmExecutor(
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
        with patch('analysis_driver.quality_control.BCLValidator._all_cycles_from_interop',
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

    @patch('analysis_driver.quality_control.bcl_validation.executor.execute', return_value=Mock(join=Mock(return_value=0)))
    @patch('analysis_driver.quality_control.BCLValidator.call_bcl_check')
    def test_check_bcls(self, mocked_check_bcls, mocked_execute):
        with patch('analysis_driver.quality_control.BCLValidator._all_cycles_from_interop',
                   return_value=[1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]):
            with patch('analysis_driver.quality_control.bcl_validation.time.sleep'):
                self.val.dataset.is_sequencing = Mock(side_effect=[True, True, False])
                self.val.check_bcls()
                assert mocked_check_bcls.called is True
                assert mocked_check_bcls.call_count == 3
