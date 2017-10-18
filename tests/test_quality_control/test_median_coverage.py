from analysis_driver.quality_control import SamtoolsDepth
from tests.test_quality_control.qc_tester import QCTester
from unittest.mock import patch


class TestSamtoolsDepth(QCTester):
    exp_cmd = ('path/to/samtools depth -a -a -q 0 -Q 0 testfile.bam | '
               'awk -F "\t" \'{array[$1"\t"$3]+=1} END{for (val in array){print val"\t"array[val]}}\' | '
               'sort -T path/to/jobs/test_sample -k 1,1 -nk 2,2 > testfile.depth')

    @patch('egcg_core.executor.execute', autospec=True)
    def test_run_samtools_depth(self, mocked_execute):
        g = SamtoolsDepth(dataset=self.dataset, bam_file='testfile.bam')
        with patch('analysis_driver.quality_control.median_coverage.find_file', new=self.fake_find_file):
            g._run()
        mocked_execute.assert_called_once_with(
            self.exp_cmd, job_name='samtoolsdepth', working_dir='path/to/jobs/test_sample', cpus=1, mem=6
        )
