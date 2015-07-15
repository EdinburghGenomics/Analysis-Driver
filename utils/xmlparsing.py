import os.path
import xml.etree.ElementTree as ET
#from util.logger import AppLogger


class RunInfo:
    def __init__(self, data_dir):
        run_info = os.path.join(data_dir, 'RunInfo.xml')
        self.root = ET.parse(run_info).getroot()
        reads = self.root.find('Run/Reads').getchildren()

        self.mask = Mask()
        for read in reads:
            self.mask.add(read)
        self.mask.validate()

        barcode_reads = self.mask.index_lengths

        if len(barcode_reads):
            self.barcode_len = int(barcode_reads[0])
        else:
            pass  # RunInfo is created with no barcode_len property


class Mask:
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
        self.reads.append(read)
        if self._is_indexed_read(read):
            assert (read.attrib == self.barcode_len or self.barcode_len is None)
            self.barcode_len = read.attrib['NumCycles']

    def validate(self):
        assert not self._is_indexed_read(self.reads[0])
        assert not self._is_indexed_read(self.reads[-1])
        for index in self.indexes:
            assert self._is_indexed_read(index), str([x.attrib for x in self.indexes])

    def tostring(self, sample_sheet_barcode_len):
        mask = 'y' + str(self._num_cycles(self.upstream_read) - 1) + 'n,'

        for i in self.index_lengths:
            diff = int(i) - sample_sheet_barcode_len
            mask += 'i' + str(sample_sheet_barcode_len) + 'n'*diff + ','

        mask += 'y' + str(self._num_cycles(self.downstream_read) - 1) + 'n'

        return mask

    def _is_indexed_read(self, read):
        if read.attrib['IsIndexedRead'] == 'Y':
            return True
        elif read.attrib['IsIndexedRead'] == 'N':
            return False
        else:
            raise ValueError('Invalid IsIndexedRead parameter: ' + read.attrib['IsIndexedRead'])

    def _num_cycles(self, read):
        return int(read.attrib['NumCycles'])



# Unit Test
if __name__ == '__main__':
    file_path = "/home/U008/lcebaman/scripts/data/RunInfo.xml"
    print(file_path)

    #mask = get_mask(file_path)
    #print(mask)
