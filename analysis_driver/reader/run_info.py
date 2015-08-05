import os.path
import xml.etree.ElementTree as eT
from analysis_driver.util import AppLogger


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
        self.root = eT.parse(run_info).getroot()
        reads = self.root.find('Run/Reads').getchildren()

        # Populate a Mask object with Read XML entities
        self.mask = Mask()
        for read in reads:
            self.mask.add(read)
        self.mask.validate()

        barcode_reads = self.mask.index_lengths

        if len(barcode_reads):
            self.debug('Barcode reads: ' + str(barcode_reads))
            self.barcode_len = int(barcode_reads[0])
        else:
            self.warn('RunInfo.xml has no barcode reads')


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

    Public properties:
        reads: All contained Read entities from RunInfo.xml
        barcode_len: The read length of the barcode to be added
        upstream_read: The read that will appear on the left hand side of the final mask
        downstream_read: Same, but for the right hand side.
        indexes: A list of all reads in between

    Public methods:
        add()
        validate()
        tostring()
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
        return self.reads[1:len(self.reads)-1]

    @property
    def index_lengths(self):
        return [read.attrib['NumCycles'] for read in self.indexes]

    def add(self, read):
        """
        Add a Read entity to self.reads. If Read is a barcode, assert that its length is consistent with the
        Reads already contained.
        :param xml.etree.ElementTree.Element read: A read entity from RunInfo.xml
        :return: None
        """
        self.reads.append(read)
        if self._is_indexed_read(read):
            assert (read.attrib == self.barcode_len or self.barcode_len is None)
            self.barcode_len = read.attrib['NumCycles']

    def validate(self):
        """
        Ensure that the first and last items of self.reads are not barcodes, and that all others are.
        :return: True if successful. Raises AssertionErrors if not.
        :rtype: bool
        :raises: AssertionError
        """
        assert not self._is_indexed_read(self.reads[0])
        assert not self._is_indexed_read(self.reads[-1])
        for index in self.indexes:
            assert self._is_indexed_read(index), str([x.attrib for x in self.indexes])
        return True

    def tostring(self, sample_sheet_barcode_len):
        """
        Translate self.reads to a formatted string that can be passed as an argument to bcl2fastq.
        :param int sample_sheet_barcode_len: An expected barcode length from SampleSheet.csv.
        :return: A bcl2fastq mask, e.g. 'y150n,i6,y150n'
        :rtype: str
        """
        mask = 'y' + str(self._num_cycles(self.upstream_read) - 1) + 'n,'

        for i in self.index_lengths:
            diff = int(i) - sample_sheet_barcode_len
            mask += 'i' + str(sample_sheet_barcode_len) + 'n'*diff + ','

        mask += 'y' + str(self._num_cycles(self.downstream_read) - 1) + 'n'

        return mask

    @staticmethod
    def _is_indexed_read(read):
        """
        Tell whether a Read is a barcode, based on its IsIndexedRead attribute.
        :param xml.etree.ElementTree.Element read: A Read entity from self.reads
        :return: True if the read's IsIndexedRead is 'Y', False if 'N'.
        :raises: ValueError if IsIndexedRead is not 'Y' or 'N'
        """
        if read.attrib['IsIndexedRead'] == 'Y':
            return True
        elif read.attrib['IsIndexedRead'] == 'N':
            return False
        else:
            raise ValueError('Invalid IsIndexedRead parameter: ' + read.attrib['IsIndexedRead'])

    @staticmethod
    def _num_cycles(read):
        """
        Translate a Read's NumCycles attrib into a read length
        :param xml.etree.ElementTree.Element read: A Read from self.reads
        :return: The Read's NumCycles attrib as a read length.
        :rtype: int
        """
        return int(read.attrib['NumCycles'])
