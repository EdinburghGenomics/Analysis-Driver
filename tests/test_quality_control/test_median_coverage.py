from analysis_driver.quality_control import GatkDepthofCoverage
from unittest.mock import patch

class TestGatkDepthOfCoverage():

    def test_get_gatk_depthofcoverage_command(self):
        bam_file = 'testfile.bam'
        sample_id = 'test_sample'
        g = GatkDepthofCoverage(sample_id = sample_id, bam_file = bam_file)
        my_gatk_depthofcoverage_command, my_gatk_depthofcoverage_outfile = g._get_gatk_depthofcoverage_command()
        assert my_gatk_depthofcoverage_command == 'java -jar GenomeAnalysisTK.jar ' \
                                                  '-T DepthOfCoverage ' \
                                                  '-R path/to/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa ' \
                                                  '-o testfile.depthofcoverage ' \
                                                  '-I testfile.bam'

        assert my_gatk_depthofcoverage_outfile == 'testfile.depthofcoverage'


    @patch('analysis_driver.executor.execute', autospec=True)
    def test_run_gatk_depthofcoverage(self, mocked_execute):
        bam_file = 'testfile.bam'
        sample_id = 'test_sample'
        g = GatkDepthofCoverage(sample_id = sample_id, bam_file = bam_file)
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        run_gatk_depthofcoverage = g._run_gatk_depthofcoverage()
        mocked_execute.assert_called_once_with(['java -jar GenomeAnalysisTK.jar '
                                                  '-T DepthOfCoverage '
                                                  '-R path/to/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa '
                                                  '-o testfile.depthofcoverage '
                                                  '-I testfile.bam'],
                                               job_name='depthofcoverage',
                                               working_dir='path/to/jobs/test_sample',
                                               cpus=2,
                                               mem=10
                                               )