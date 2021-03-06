import os
import xml.etree.ElementTree as eT
from analysis_driver.reader.run_info import RunInfo, Reads
from tests.test_analysisdriver import TestAnalysisDriver


def new_read(number, num_cycles, is_indexed_read):
    return eT.Element(
        'Read', attrib={'Number': str(number), 'NumCycles': str(num_cycles), 'IsIndexedRead': is_indexed_read}
    )


def elements_equal(observed, expected):
    return sorted(observed.items()) == sorted(expected.items())


class TestRunInfo(TestAnalysisDriver):
    def setUp(self):
        self.run_info = RunInfo(self.assets_path)

    def test_flowcellname(self):
        assert self.run_info.flowcell_name == 'H5VJMCCXX'

    def test_tiles(self):
        assert self.run_info.tiles == [
            '1_1101', '2_1101', '3_1101', '4_1101', '5_1101', '6_1101', '7_1101', '8_1101'
        ]


class TestMask(TestAnalysisDriver):
    def setUp(self):
        barcoded_assets_path = os.path.join(self.assets_path, 'test_runs', 'barcoded_run')
        barcodeless_assets_path = os.path.join(self.assets_path, 'test_runs', 'barcodeless_run')

        self.barcoded_run_info = RunInfo(barcoded_assets_path)
        self.barcodeless_run_info = RunInfo(barcodeless_assets_path)

    def test_add_reads(self):
        reads = [new_read(1337, 75, 'N'), new_read(1338, 15, 'Y'), new_read(1339, 75, 'N')]
        r = Reads(reads)
        assert reads == r.reads
        assert r.barcode_len == 15

    def test_upstream_read(self):
        assert elements_equal(self.barcoded_run_info.reads.upstream_read, new_read(1, 151, 'N'))
        assert elements_equal(self.barcodeless_run_info.reads.upstream_read, new_read(1, 151, 'N'))

    def test_downstream_read(self):
        assert elements_equal(self.barcoded_run_info.reads.upstream_read, new_read(1, 151, 'N'))
        assert elements_equal(self.barcodeless_run_info.reads.upstream_read, new_read(1, 151, 'N'))

    def test_indexes(self):
        assert elements_equal(self.barcoded_run_info.reads.indexes[0], new_read(2, 8, 'Y'))
        assert self.barcodeless_run_info.reads.indexes == []

    def test_index_lengths(self):
        assert self.barcoded_run_info.reads.index_lengths == [8]
        assert self.barcodeless_run_info.reads.index_lengths == []

    def test_generate_mask(self):
        assert self.barcoded_run_info.reads.generate_mask(8) == 'y150n,i8,y150n'
        assert self.barcoded_run_info.reads.generate_mask(6) == 'y150n,i6nn,y150n'
        assert self.barcoded_run_info.reads.generate_mask(0) == 'y150n,nnnnnnnn,y150n'
        assert self.barcodeless_run_info.reads.generate_mask(8) == 'y150n,y150n'
        assert self.barcodeless_run_info.reads.generate_mask(0) == 'y150n,y150n'
