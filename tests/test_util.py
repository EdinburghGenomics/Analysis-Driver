from unittest.mock import patch, Mock
from analysis_driver.util import find_all_fastq_pairs_for_lane, get_ranges, get_trim_values_for_bad_cycles


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
