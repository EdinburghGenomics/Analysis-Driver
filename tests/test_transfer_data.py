import os
import shutil
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import transfer_data
from analysis_driver.config import OutputFileConfiguration


class TestTransferData(TestAnalysisDriver):
    data_output = os.path.join(TestAnalysisDriver.assets_path, 'data_output')
    sample_id = '10015AT0001'

    def setUp(self):
        self.link_dir = os.path.join(self.data_output, 'linked_output_files')
        os.makedirs(self.link_dir, exist_ok=True)
        self.output_cfg = OutputFileConfiguration('non_human_qc')

    def tearDown(self):
        if os.path.isdir(self.link_dir):
            shutil.rmtree(self.link_dir)

    def test_create_links(self):
        list_of_linked_files = transfer_data.create_output_links(
            self.data_output, self.output_cfg, self.link_dir, sample_id=self.sample_id, user_sample_id=self.sample_id
        )

        output_files = os.path.join(self.data_output, 'linked_output_files')

        expected_outputs = ['10015AT0001.depth', '10015AT0001_R1_fastqc.html', '10015AT0001_R1_fastqc.zip',
                            '10015AT0001_R1_screen.txt', '10015AT0001_R2_fastqc.html',
                            '10015AT0001_R2_fastqc.zip', '10015AT0001_filter_snp.vcf.stats',
                            'samtools_stats.txt', 'taxa_identified.json']
        assert sorted(os.listdir(output_files)) == expected_outputs == sorted(
            os.path.basename(f) for f in list_of_linked_files
        )
