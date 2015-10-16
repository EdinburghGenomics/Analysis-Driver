__author__ = 'mwham'
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.reader.run_info import RunInfo, Mask
import xml.etree.ElementTree as eT


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
        self.run_info_helper = RunInfo(self.assets_path)
        self.mask_helper = self.run_info_helper.mask

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
        assert self.mask_helper.upstream_read.attrib['NumCycles'] == '151'

    def test_downstream_read(self):
        assert self.mask_helper.upstream_read.attrib['NumCycles'] == '151'

    def test_indexes(self):
        assert [x.attrib['NumCycles'] for x in self.run_info_helper.mask.indexes] == ['6']

    def test_index_lengths(self):
        assert self.mask_helper.index_lengths == [6]

    def test_validate(self):
        assert self.mask_helper.validate()


