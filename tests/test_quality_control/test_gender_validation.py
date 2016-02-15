from unittest.mock import patch
from analysis_driver.quality_control.gender_validation import GenderValidation
from tests.test_analysisdriver import TestAnalysisDriver

__author__ = 'tcezard'


class TestGenderValidation(TestAnalysisDriver):

    @patch('analysis_driver.executor.execute')
    def test__gender_call(self, mocked_execute):
        validator = GenderValidation(sample_id='test_sample', vcf_file='path/to/test/vcf')
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        return_code = validator._gender_call()
        assert return_code == 0
        assert mocked_execute.call_args[0][0][0].startswith('cat')
        instance.join.return_value = 1
        return_code = validator._gender_call()
        assert return_code == 1
        validator = GenderValidation(sample_id='test_sample', vcf_file='path/to/test/vcf.gz')
        return_code = validator._gender_call()
        assert mocked_execute.call_args[0][0][0].startswith('zcat')




