import glob
import os
import shutil
from unittest.mock import patch, Mock
from analysis_driver.util import find_all_fastq_pairs_for_lane, get_ranges, get_trim_values_for_bad_cycles
from analysis_driver.util.helper_functions import prepend_path_to_data_files, split_in_chunks, merge_lane_directories
from tests import TestAnalysisDriver


def test_get_range():
    list_int = [1, 2, 3, 4, 6, 7]
    assert list(get_ranges(list_int)) == [(1, 4), (6, 7)]


def test_find_fastq_pairs():
    fqs = ('test_L001_R1.fastq.gz', 'test_L001_R2.fastq.gz', 'test_L002_R1.fastq.gz', 'test_L002_R2.fastq.gz')
    fake_walk = (
        ('run_dir', ['a_subdir', 'another_subdir'], fqs),
        ('run_dir/a_subdir', [], fqs),
        ('run_dir/another_subdir', [], fqs[:2])
    )
    with patch('analysis_driver.util.os.walk', return_value=fake_walk):
        assert find_all_fastq_pairs_for_lane('a_location', 1) == [
            ('run_dir/a_subdir/test_L001_R1.fastq.gz', 'run_dir/a_subdir/test_L001_R2.fastq.gz'),
            ('run_dir/another_subdir/test_L001_R1.fastq.gz', 'run_dir/another_subdir/test_L001_R2.fastq.gz'),
            ('run_dir/test_L001_R1.fastq.gz', 'run_dir/test_L001_R2.fastq.gz')
        ]
        assert find_all_fastq_pairs_for_lane('a_location', 2) == [
            ('run_dir/a_subdir/test_L002_R1.fastq.gz', 'run_dir/a_subdir/test_L002_R2.fastq.gz'),
            ('run_dir/test_L002_R1.fastq.gz', 'run_dir/test_L002_R2.fastq.gz')
        ]


def test_trim_values_for_bad_cycles():
    run_info = Mock(reads=Mock(
        upstream_read=Mock(attrib={'NumCycles': '151'}),
        downstream_read=Mock(attrib={'NumCycles': '151'}),
        index_lengths=[8]
    ))
    bad_cycle_list = [310, 308, 307, 309, 101]
    assert get_trim_values_for_bad_cycles(bad_cycle_list, run_info) == (None, 147)

    bad_cycle_list = [308, 307, 309, 151, 150]
    assert get_trim_values_for_bad_cycles(bad_cycle_list, run_info) == (149, None)

    bad_cycle_list = [310, 309]
    assert get_trim_values_for_bad_cycles(bad_cycle_list, run_info) == (None, 149)

    bad_cycle_list = None
    assert get_trim_values_for_bad_cycles(bad_cycle_list, run_info) == (None, None)

    run_info = Mock(reads=Mock(
        upstream_read=Mock(attrib={'NumCycles': '151'}),
        downstream_read=Mock(attrib={'NumCycles': '151'}),
        index_lengths=[]
    ))

    bad_cycle_list = [302, 301]
    assert get_trim_values_for_bad_cycles(bad_cycle_list, run_info) == (None, 149)


class TestHelperFunctions(TestAnalysisDriver):

    def test_prepend_path_to_data_files(self):
        file_paths = {'f1': 'path/file.1', 'file_set': {'f2': 'path/file.2', 'f3': 'path/file.3'}}
        full_file_paths = {
            'f1': '/full/path/file.1',
            'file_set': {'f2': '/full/path/file.2', 'f3': '/full/path/file.3'}
        }
        # Keeps relative path when appending nothing
        assert prepend_path_to_data_files('', file_paths) == file_paths

        # Prepend path if provided
        assert prepend_path_to_data_files('/full', file_paths) == full_file_paths

        # Does not prepend anything on already full path
        assert prepend_path_to_data_files('/fuller', full_file_paths) == full_file_paths

    def test_split_in_chunks(self):
        chunks = split_in_chunks(total_length=15, chunksize=20, zero_based=True, end_inclusive=True)
        assert chunks == [(0, 14)]

        chunks = split_in_chunks(total_length=15, chunksize=20, zero_based=False, end_inclusive=True)
        assert chunks == [(1, 15)]

        chunks = split_in_chunks(total_length=15, chunksize=20, zero_based=True, end_inclusive=False)
        assert chunks == [(0, 15)]

        chunks = split_in_chunks(total_length=15, chunksize=20, zero_based=False, end_inclusive=False)
        assert chunks == [(1, 16)]

        chunks = split_in_chunks(total_length=122, chunksize=20, zero_based=True, end_inclusive=True)
        assert chunks == [(0, 19), (20, 39), (40, 59), (60, 79), (80, 99), (100, 119), (120, 121)]

        chunks = split_in_chunks(total_length=122, chunksize=20, zero_based=False, end_inclusive=True)
        assert chunks == [(1, 20), (21, 40), (41, 60), (61, 80), (81, 100), (101, 120), (121, 122)]

        chunks = split_in_chunks(total_length=122, chunksize=20, zero_based=True, end_inclusive=False)
        assert chunks == [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100), (100, 120), (120, 122)]

        chunks = split_in_chunks(total_length=122, chunksize=20, zero_based=False, end_inclusive=False)
        assert chunks == [(1, 21), (21, 41), (41, 61), (61, 81), (81, 101), (101, 121), (121, 123)]

    def test_merge_lane_directories(self):
        fastq_dir = os.path.join(self.assets_path, 'fastq_split')
        os.makedirs(fastq_dir, exist_ok=True)
        run_elements = [
            {'project_id': 'a_project', 'lane': '1', 'sample_id': 'sample1'},
            {'project_id': 'a_project', 'lane': '1', 'sample_id': 'sample2'},
            {'project_id': 'a_project', 'lane': '2', 'sample_id': 'sample1'},
            {'project_id': 'a_project', 'lane': '2', 'sample_id': 'sample2'}
        ]
        for re in run_elements:
            sample_dir = os.path.join(fastq_dir, 'lane_%s' % re['lane'], re['project_id'], re['sample_id'])
            os.makedirs(sample_dir, exist_ok=True)
            self._touch(os.path.join(sample_dir, '%s_L00%s_R1_001.fastq.gz' % (re['sample_id'], re['lane'])))
            self._touch(os.path.join(sample_dir, '%s_L00%s_R2_001.fastq.gz' % (re['sample_id'], re['lane'])))

        # Add unassigned
        self._touch(os.path.join(fastq_dir, 'lane_1', 'Undetermined_S0_L001_R1_001.fastq.gz'))
        self._touch(os.path.join(fastq_dir, 'lane_1', 'Undetermined_S0_L001_R2_001.fastq.gz'))

        # Add stats files
        os.makedirs(os.path.join(fastq_dir, 'lane_1', 'Stats'), exist_ok=True)
        os.makedirs(os.path.join(fastq_dir, 'lane_2', 'Stats'), exist_ok=True)
        self._touch(os.path.join(fastq_dir, 'lane_1', 'Stats', 'Stats.json'))
        self._touch(os.path.join(fastq_dir, 'lane_2', 'Stats', 'Stats.json'))

        merge_lane_directories(fastq_dir, run_elements)

        #shutil.rmtree(fastq_dir)
