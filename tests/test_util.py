from unittest.mock import patch, Mock

import os

from analysis_driver.util import find_all_fastq_pairs_for_lane, get_ranges, get_trim_values_for_bad_cycles, \
    create_regions_of_genome


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


def test_create_regions_of_genome():
    assets_path = os.path.join(os.path.dirname(__file__), 'assets')
    reference_genome = os.path.join(assets_path, 'genomes', 'genome.fa')
    bed_file = os.path.join(assets_path, 'genomes', 'genome_region.bed')
    expected_bed_lines = ['ref4	0	12500\n', 'ref5	0	8500\n', 'ref7	0	6800\n', 'ref8	0	3800\n',
                          'ref2	0	3600\n', 'ref1	0	2500\n', 'ref6	0	1800\n', 'ref3	0	1100\n']
    create_regions_of_genome(reference_genome, output_bed=bed_file)
    with open(bed_file) as open_bed:
        assert open_bed.readlines() == expected_bed_lines

    create_regions_of_genome(reference_genome, output_bed=bed_file, max_nb_contig=4)
    with open(bed_file) as open_bed:
        assert open_bed.readlines() == expected_bed_lines[:4]

    create_regions_of_genome(reference_genome, output_bed=bed_file, regions_size=35000)
    with open(bed_file) as open_bed:
        assert open_bed.readlines() == expected_bed_lines[:5]
