import os
import xml.etree.ElementTree as eT
from analysis_driver.reader.run_info import RunInfo, Mask
from tests.test_analysisdriver import TestAnalysisDriver


def new_read(number, num_cycles, is_indexed_read):
    return eT.Element(
        'Read', attrib={'Number': str(number), 'NumCycles': str(num_cycles), 'IsIndexedRead': is_indexed_read}
    )


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

        self.barcoded_run_info_helper = RunInfo(barcoded_assets_path)
        self.barcoded_mask_helper = self.barcoded_run_info_helper.mask

        self.barcodeless_run_info_helper = RunInfo(barcodeless_assets_path)
        self.barcodeless_mask_helper = self.barcodeless_run_info_helper.mask

    def test_add(self):
        new_element = new_read(1337, 75, 'N')
        self.mask.add(new_element)
        assert new_element in self.mask.reads
        assert self.mask.barcode_len is None

        new_element = new_read(1338, 15, 'Y')
        self.mask.add(new_element)
        assert new_element in self.mask.reads
        assert self.mask.barcode_len == 15

        new_element = new_read(1339, 75, 'N')
        self.mask.add(new_element)
        assert new_element in self.mask.reads
        assert self.mask.barcode_len == 15



    def test_upstream_read(self):
        assert self.barcoded_mask_helper.upstream_read.attrib['NumCycles'] == '151'
        assert self.barcodeless_mask_helper.upstream_read.attrib['NumCycles'] == '151'

    def test_downstream_read(self):
        assert self.barcoded_mask_helper.upstream_read.attrib['NumCycles'] == '151'
        assert self.barcodeless_mask_helper.upstream_read.attrib['NumCycles'] == '151'

    def test_indexes(self):
        assert [x.attrib['NumCycles'] for x in self.barcoded_run_info_helper.mask.indexes] == ['6']
        #self.assertRaises(ValueError, self.barcodeless_run_info_helper.mask.indexes)

    def test_index_lengths(self):
        assert self.barcoded_mask_helper.index_lengths == [6]
    #    self.assertRaises(AnalysisDriverError, self.barcodeless_mask_helper.index_lengths)

    #def test_validate_barcoded(self):
    #    assert self.barcoded_mask_helper.validate_barcoded()

    #def test_validate_barcodeless(self):
    #    assert self.barcodeless_mask_helper.validate_barcodeless()