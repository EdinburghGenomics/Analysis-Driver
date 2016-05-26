from unittest.mock import patch
from analysis_driver.quality_control.gender_validation import GenderValidation
from tests.test_quality_control.qc_tester import QCTester


class TestGenderValidation(QCTester):
    @patch('analysis_driver.dataset.rest_communication')
    @patch('analysis_driver.executor.execute')
    def test__gender_call(self, mocked_execute, mocked_rest):
        validator = GenderValidation(self.dataset, working_dir='test_sample', vcf_file='path/to/test/vcf')
        instance = mocked_execute.return_value

        instance.join.return_value = 0
        return_code = validator._gender_call()
        assert return_code == 0
        command = mocked_execute.call_args[0][0]
        assert command.startswith('cat')
        assert len(command.split(' | ')) == 5

        instance.join.return_value = 1
        return_code = validator._gender_call()
        assert return_code == 1

        validator = GenderValidation(self.dataset, working_dir='test_sample', vcf_file='path/to/test/vcf.gz')
        validator._gender_call()
        assert mocked_execute.call_args[0][0].startswith('zcat')
