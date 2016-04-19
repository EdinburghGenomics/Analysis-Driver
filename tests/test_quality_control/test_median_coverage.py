from analysis_driver.quality_control import SamtoolsDepth
from tests.test_quality_control.qc_tester import QCTester
from unittest.mock import patch

class TestSamtoolsDepth(QCTester):

    def test_get_samtools_depth_command(self):
        bam_file = 'testfile.bam'
        working_dir = 'test_sample'
        g = SamtoolsDepth(self.dataset, bam_file = bam_file, working_dir = working_dir)
        my_samtools_depth_command, my_samtools_depth_outfile = g._get_samtools_depth_command()
        assert my_samtools_depth_command == """path/to/samtools depth -a -a -q 0 -Q 0 testfile.bam | awk -F "\t" '{array[$1"\t"$3]+=1} END{for (val in array){print val"\t"array[val]}}' | sort -k 1,1 -nk 2,2 > test_sample/samtoolsdepth.txt"""

        assert my_samtools_depth_outfile == 'test_sample/samtoolsdepth.txt'


    @patch('analysis_driver.executor.execute', autospec=True)
    def test_run_samtools_depth(self, mocked_execute):
        bam_file = 'testfile.bam'
        working_dir = 'test_sample'
        g = SamtoolsDepth(self.dataset, bam_file = bam_file, working_dir = working_dir)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_samtools_depth = g._run_samtools_depth()
        mocked_execute.assert_called_once_with(["""path/to/samtools depth -a -a -q 0 -Q 0 testfile.bam | awk -F "\t" '{array[$1"\t"$3]+=1} END{for (val in array){print val"\t"array[val]}}' | sort -k 1,1 -nk 2,2 > test_sample/samtoolsdepth.txt"""],
                                               job_name='samtoolsdepth',
                                               working_dir='test_sample',
                                               cpus=1,
                                               mem=6
                                               )


