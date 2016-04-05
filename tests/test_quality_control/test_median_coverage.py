from analysis_driver.quality_control import SamtoolsDepth
from tests.test_quality_control.qc_tester import QCTester
from unittest.mock import patch

class TestSamtoolsDepth(QCTester):

    def test_get_samtools_depth_command(self):
        bam_file = 'testfile.bam'
        working_dir = 'test_sample'
        g = SamtoolsDepth(self.dataset, bam_file = bam_file, working_dir = working_dir)
        my_samtools_depth_command, my_samtools_depth_outfile = g._get_samtools_depth_command()
        assert my_samtools_depth_command == 'path/to/samtools depth -a -a -q 0 -Q 0 testfile.bam > testfile.depth'

        assert my_samtools_depth_outfile == 'testfile.depth'


    @patch('analysis_driver.executor.execute', autospec=True)
    def test_run_samtools_depth(self, mocked_execute):
        bam_file = 'testfile.bam'
        working_dir = 'test_sample'
        g = SamtoolsDepth(self.dataset, bam_file = bam_file, working_dir = working_dir)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_samtools_depth = g._run_samtools_depth()
        mocked_execute.assert_called_once_with(['path/to/samtools depth -a -a -q 0 -Q 0 testfile.bam > testfile.depth'],
                                               job_name='samtoolsdepth',
                                               working_dir='test_sample',
                                               cpus=2,
                                               mem=10
                                               )

    @patch('analysis_driver.quality_control.SamtoolsDepth._run_samtools_depth')
    def test_get_depth_histogram_command(self, mocked_depthfile):
        bam_file = 'testfile.bam'
        working_dir = 'test_sample'
        g = SamtoolsDepth(self.dataset, bam_file = bam_file, working_dir = working_dir)
        mocked_depthfile.return_value = 'test_sample.depth'
        depth_histogram_command, depth_histogram_outfile = g._get_depth_histogram_command()
        assert depth_histogram_command == "awk -F '\t' '{print $3}' test_sample.depth | sort | uniq -c | sort -nr > test_sample.hist"


