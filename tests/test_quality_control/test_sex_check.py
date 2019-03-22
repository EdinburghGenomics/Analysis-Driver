from unittest.mock import patch
from analysis_driver.quality_control.sex_check import SexCheck
from tests.test_quality_control.qc_tester import QCTester


class TestSexCheck(QCTester):
    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        validator = SexCheck(dataset=self.dataset, vcf_file='path/to/test/vcf')

        with patch('analysis_driver.quality_control.sex_check.util.find_file', new=self.fake_find_file):
            validator._run()

            command = mocked_execute.call_args[0][0]
            assert command.startswith('cat')
            assert len(command.split(' | ')) == 5

            validator = SexCheck(dataset=self.dataset, vcf_file='path/to/test/vcf.gz')
            validator._run()
            assert mocked_execute.call_args[0][0].startswith('zcat')
