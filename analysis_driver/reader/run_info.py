import os.path
from xml.etree import ElementTree
from analysis_driver.exceptions import AnalysisDriverError
from egcg_core.app_logging import AppLogger


class RunInfo(AppLogger):
    """
    Represents a RunInfo.xml file. It uses xml.etree to read the file from disk and populate a Mask
    object.
    """

    def __init__(self, data_dir):
        """
        :param str data_dir: A file path to the input_data folder containing RunInfo.xml
        """
        run_info = os.path.join(data_dir, 'RunInfo.xml')
        root = ElementTree.parse(run_info).getroot()
        reads = root.find('Run/Reads').getchildren()

        self.reads = Reads(reads)
        self.flowcell_name = root.find('Run/Flowcell').text
        self.tiles = [e.text for e in root.find('Run/FlowcellLayout/TileSet/Tiles')]


class Reads:
    """
    Represents a series of Read entities from RunInfo.xml (contained in self.reads) and stores the length of
    all the barcodes (indexed Reads), ensuring they are the same length. Reads are stored in order. The
    NumCycles and IsIndexedRead attributes are translated into a bcl2fastq mask.
    """

    def __init__(self, reads):
        self.reads = []
        self.barcode_len = None
        self.indexes = []

        self._add_reads(reads)
        self.upstream_read = self.reads[0]
        self.downstream_read = self.reads[-1]
        self.index_lengths = [self.num_cycles(i) for i in self.indexes]

    @property
    def has_barcodes(self):
        return self.barcode_len is not None

    def _add_reads(self, reads):
        """
        Add a Read entity to self.reads and if it is a barcode, assert that its length is consistent with the
        Reads already contained.
        :param et.Element reads: Read entities from RunInfo.xml
        """
        for r in reads:
            self.reads.append(r)
            if self._is_indexed_read(r):
                self.indexes.append(r)
                assert r.attrib['NumCycles'] == self.barcode_len or self.barcode_len is None
                self.barcode_len = int(r.attrib['NumCycles'])

        if not self.reads:
            raise AnalysisDriverError('No reads found in RunInfo.xml')

    @staticmethod
    def num_cycles(read):
        return int(read.attrib['NumCycles'])

    @staticmethod
    def _is_indexed_read(read):
        """Translate IsIndexedRead from "Y"/"N" to True/False."""
        if read.attrib['IsIndexedRead'] not in 'YN':
            raise AnalysisDriverError('Invalid IsIndexedRead parameter: ' + read.attrib['IsIndexedRead'])
        return read.attrib['IsIndexedRead'] == 'Y'

    def generate_mask(self, samples_barcode_len):
        """
        Translate:
            <Read IsIndexedRead=N Number=1 NumCycles=151/>
            <Read IsIndexedRead=Y Number=2 NumCycles=8/>
            <Read IsIndexedRead=N Number=3 NumCycles=151/>
        to 'y150n,i8,y150n'. If the sample sheet says the barcode is shorter, trailing 'n's are added, e.g.
        'y150n,i6nn,y150n'.
        """
        out = ['y' + str(self.num_cycles(self.upstream_read) - 1) + 'n']

        for i in self.index_lengths:
            diff = i - samples_barcode_len
            out.append('i' + str(samples_barcode_len) + 'n' * diff)

        out.append('y' + str(self.num_cycles(self.downstream_read) - 1) + 'n')
        return ','.join(out)
