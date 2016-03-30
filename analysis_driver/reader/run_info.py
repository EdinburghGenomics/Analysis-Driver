import os.path
from xml.etree import ElementTree
from analysis_driver.app_logging import AppLogger


class RunInfo(AppLogger):
    """
    Represents a RunInfo.xml file. It uses xml.etree to read the file from disk and populate a Mask
    object.

    Public properties:
        mask: A Mask object
        barcode_len: The shared read length for all barcode reads.
    """

    def __init__(self, data_dir):
        """
        :param str data_dir: A file path to the input_data folder containing RunInfo.xml
        """
        run_info = os.path.join(data_dir, 'RunInfo.xml')
        root = ElementTree.parse(run_info).getroot()
        reads = root.find('Run/Reads').getchildren()

        # Populate a Mask object with Read XML entities
        self.mask = Mask()
        for read in reads:
            self.debug('Adding read: ' + str(read.attrib))
            self.mask.add(read)

        self.flowcell_name = root.find('Run/Flowcell').text


class Mask:
    """
    Represents a series of Read XML entities from RunInfo.xml, contained in self.reads. Also stores the
    length of the barcode Reads and ensures that they are all the same.

    RunInfo.xml contains the following:

    <Reads>
      <Read Number="1" NumCycles="151" IsIndexedRead="N" />
      <Read Number="2" NumCycles="6" IsIndexedRead="Y" />
      <Read Number="3" NumCycles="151" IsIndexedRead="N" />
    </Reads>

    The Reads are stored in order. NumCycles represents the read length and IsIndexedRead represents
    whether the Read is a barcode. These are translated to a string to be passed to bcl2fastq as a mask. The
    example here would be translated to a mask of: 'y150n,i6,y150n'.
    """

    def __init__(self):
        self.reads = []
        self.barcode_len = None

    @property
    def upstream_read(self):
        return self.reads[0]

    @property
    def downstream_read(self):
        return self.reads[-1]

    @property
    def indexes(self):
        return [r for r in self.reads if self._is_indexed_read(r)]

    @property
    def index_lengths(self):
        return [self.num_cycles(read) for read in self.indexes]

    @property
    def has_barcodes(self):
        return self.barcode_len is not None

    def add(self, read):
        """
        Add a Read entity to self.reads and if it is a barcode, assert that its length is consistent with the
        Reads already contained.
        :param et.Element read: A read entity from RunInfo.xml
        """
        self.reads.append(read)
        if self._is_indexed_read(read):
            assert (read.attrib == self.barcode_len or self.barcode_len is None)
            self.barcode_len = int(read.attrib['NumCycles'])

    @staticmethod
    def num_cycles(read):
        """
        Return a RunInfo.xml Read's NumCycles attrib as a read length
        """
        return int(read.attrib['NumCycles'])

    @staticmethod
    def _is_indexed_read(read):
        """
        Tell whether a RunInfo.xml Read is a barcode, based on its IsIndexedRead attribute.
        """
        if read.attrib['IsIndexedRead'] == 'Y':
            return True
        elif read.attrib['IsIndexedRead'] == 'N':
            return False
        else:
            raise ValueError('Invalid IsIndexedRead parameter: ' + read.attrib['IsIndexedRead'])
