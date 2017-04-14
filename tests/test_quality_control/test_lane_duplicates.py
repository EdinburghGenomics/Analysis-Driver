from unittest.mock import patch
from analysis_driver.quality_control.lane_duplicates import WellDuplicates
from tests.test_quality_control.qc_tester import QCTester


class TestWellDuplicates(QCTester):
    @patch('egcg_core.executor.execute')
    def test_well_duplicates(self, mocked_execute):
        welldup = WellDuplicates(
            dataset=self.run_dataset,
            output_directory='test_run/fastq',
            run_directory='path/to/run'
        )
        welldup._run()
        mocked_execute.assert_called_once_with(
            ('path/to/well_duplicate -f path/to/coord_file -r path/to/run -s hiseq_x > '
             'test_run/fastq/test_run.wellduplicate 2> test_run/fastq/test_run.wellduplicate.err'),
            cpus=1,
            job_name='welldup',
            mem=2,
            log_commands=False,
            working_dir='path/to/jobs/test_run'
        )
