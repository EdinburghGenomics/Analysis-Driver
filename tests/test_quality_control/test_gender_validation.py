from unittest.mock import patch
from analysis_driver.quality_control.gender_validation import GenderValidation
from tests.test_quality_control.qc_tester import QCTester


class TestGenderValidation(QCTester):
    @patch('egcg_core.executor.execute')
    def test_run(self, mocked_execute):
        validator = GenderValidation(dataset=self.dataset, vcf_file='path/to/test/vcf')

        validator._run()
        command = mocked_execute.call_args[0][0]
        assert command.startswith('cat')
        assert len(command.split(' | ')) == 5

        validator = GenderValidation(dataset=self.dataset, vcf_file='path/to/test/vcf.gz')
        validator._run()
        assert mocked_execute.call_args[0][0].startswith('zcat')
