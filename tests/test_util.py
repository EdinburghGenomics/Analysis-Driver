import os
import shutil
from unittest.mock import patch, Mock

from analysis_driver.util import get_ranges, convert_bad_cycle_in_trim
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver import transfer_data, util
from analysis_driver.config import OutputFileConfiguration


def patched_find_project_from_sample(sample_id):
    return patch('egcg_core.clarity.find_project_name_from_sample', return_value='proj_' + sample_id)


class TestTransferData(TestAnalysisDriver):
    data_output = os.path.join(TestAnalysisDriver.assets_path, 'data_output')
    sample_id = '10015AT0001'

    def setUp(self):
        self.link_dir = os.path.join(self.data_output, 'linked_output_files')
        os.makedirs(self.link_dir, exist_ok=True)

        self._to_dir = os.path.join(self.data_output, 'to', '')
        os.makedirs(self._to_dir, exist_ok=True)
        self._pseudo_links = os.path.join(self.data_output, 'pseudo_links')
        self.output_cfg = OutputFileConfiguration('non_human_qc')
        self._create_pseudo_links()

    def tearDown(self):
        shutil.rmtree(self._to_dir)
        shutil.rmtree(self._pseudo_links)
        if os.path.isdir(self.link_dir):
            shutil.rmtree(self.link_dir)

    def _create_pseudo_links(self):
        os.makedirs(self._pseudo_links, exist_ok=True)
        for k in self.output_cfg.content:
            f = os.path.join(
                self._pseudo_links,
                self.output_cfg.output_dir_file(k).format(
                    sample_id=self.sample_id, user_sample_id=self.sample_id
                )
            )
            open(f, 'a').close()

    def test_create_links(self):
        with patch('egcg_core.clarity.get_user_sample_name', return_value=self.sample_id):
            list_of_linked_files = transfer_data.create_links_from_bcbio(
                self.sample_id, self.data_output, self.output_cfg, self.link_dir
            )

        output_files = os.path.join(self.data_output, 'linked_output_files')

        expected_outputs = ['10015AT0001.depth', '10015AT0001_R1_fastqc.html', '10015AT0001_R1_fastqc.zip',
                            '10015AT0001_R1_screen.txt', '10015AT0001_R2_fastqc.html',
                            '10015AT0001_R2_fastqc.zip', 'samtools_stats.txt', 'taxa_identified.json']
        assert sorted(os.listdir(output_files)) == expected_outputs == sorted(
            os.path.basename(f) for f in list_of_linked_files
        )

    def test_output_sample_data(self):
        with patched_find_project_from_sample(self.sample_id), \
                patch('analysis_driver.transfer_data.archive_management.archive_directory', return_value=True):
            exit_status = transfer_data.output_sample_data(
                sample_id=self.sample_id,
                source_dir=self._pseudo_links,
                output_dir=self._to_dir
            )
        output_files = os.path.join(self._to_dir, 'proj_' + self.sample_id, self.sample_id)

        expected_outputs = ['10015AT0001.depth', '10015AT0001_R1_fastqc.html', '10015AT0001_R1_fastqc.zip',
                            '10015AT0001_R1_screen.txt', '10015AT0001_R2_fastqc.html',
                            '10015AT0001_R2_fastqc.zip', 'samtools_stats.txt', 'taxa_identified.json']

        o = sorted(os.listdir(output_files))
        assert exit_status == 0
        assert o == expected_outputs

    def test_prep_samples_cmd(self):
        cmd = util.bcbio_prepare_samples_cmd(
            self.assets_path,
            'a_sample_id',
            ['test_R1.fastq', 'test_R2.fastq'],
            user_sample_id='a_user_sample_id'
        )

        bcbio_csv = os.path.join(self.assets_path, 'samples_a_sample_id.csv')
        with open(bcbio_csv) as f:
            content = f.read()
            print(content)
            assert content == (
                'samplename,description\n'
                'test_R1.fastq,a_user_sample_id\n'
                'test_R2.fastq,a_user_sample_id\n'
            )
        os.remove(bcbio_csv)
        assert not os.path.isfile(bcbio_csv)
        expected = (
            'bcbio_prepare_samples.py --out ' +
            os.path.join(self.assets_path, 'merged') +
            ' --csv ' +
            bcbio_csv
        )
        assert expected in cmd


class TestUtils(TestAnalysisDriver):

    def test_get_range(self):

        list_int = [1, 2, 3, 4, 6, 7]
        assert list(get_ranges(list_int)) == [(1, 4), (6, 7)]


    def test_convert_bad_cycle_in_trim(self):
        run_info = Mock(reads=Mock(
            upstream_read=Mock(attrib={'NumCycles': '151'}),
            downstream_read=Mock(attrib={'NumCycles': '151'}),
            index_lengths=[8]
        ))
        bad_cycle_list = [310, 308, 307, 309, 101]
        assert convert_bad_cycle_in_trim(bad_cycle_list, run_info) == (None, 147)

        bad_cycle_list = [308, 307, 309, 151, 150]
        assert convert_bad_cycle_in_trim(bad_cycle_list, run_info) == (149, None)

        bad_cycle_list = [310, 309]
        assert convert_bad_cycle_in_trim(bad_cycle_list, run_info) == (None, 149)

        bad_cycle_list = None
        assert convert_bad_cycle_in_trim(bad_cycle_list, run_info) == (None, None)

        run_info = Mock(reads=Mock(
            upstream_read=Mock(attrib={'NumCycles': '151'}),
            downstream_read=Mock(attrib={'NumCycles': '151'}),
            index_lengths=[]
        ))

        bad_cycle_list = [302, 301]
        assert convert_bad_cycle_in_trim(bad_cycle_list, run_info) == (None, 149)


