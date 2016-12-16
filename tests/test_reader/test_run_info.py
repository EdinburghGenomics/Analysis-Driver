import os
import xml.etree.ElementTree as eT
from analysis_driver.reader.run_info import RunInfo, Mask
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


class TestMask(TestAnalysisDriver):
    def setUp(self):
        self.mask = Mask()

        barcoded_assets_path = os.path.join(self.assets_path, 'test_runs', 'barcoded_run')
        barcodeless_assets_path = os.path.join(self.assets_path, 'test_runs', 'barcodeless_run')

        self.barcoded_run_info = RunInfo(barcoded_assets_path)
        self.barcodeless_run_info = RunInfo(barcodeless_assets_path)

    def test_add(self):
        r = new_read(1337, 75, 'N')
        self.mask.add(r)
        assert r in self.mask.reads
        assert self.mask.barcode_len is None

        r = new_read(1338, 15, 'Y')
        self.mask.add(r)
        assert r in self.mask.reads
        assert self.mask.barcode_len == 15

        r = new_read(1339, 75, 'N')
        self.mask.add(r)
        assert r in self.mask.reads
        assert self.mask.barcode_len == 15

    def test_upstream_read(self):
        assert elements_equal(self.barcoded_run_info.mask.upstream_read, new_read(1, 151, 'N'))
        assert elements_equal(self.barcodeless_run_info.mask.upstream_read, new_read(1, 151, 'N'))

    def test_downstream_read(self):
        assert elements_equal(self.barcoded_run_info.mask.upstream_read, new_read(1, 151, 'N'))
        assert elements_equal(self.barcodeless_run_info.mask.upstream_read, new_read(1, 151, 'N'))

    def test_indexes(self):
        assert elements_equal(self.barcoded_run_info.mask.indexes[0], new_read(2, 8, 'Y'))
        assert self.barcodeless_run_info.mask.indexes == []

    def test_index_lengths(self):
        assert self.barcoded_run_info.mask.index_lengths == [8]
        assert self.barcodeless_run_info.mask.index_lengths == []

    # def test_validate_barcoded(self):
    #     assert self.barcoded_run_info.mask.validate_barcoded()
    #
    # def test_validate_barcodeless(self):
    #     assert self.barcodeless_run_info.mask.validate_barcodeless()
