import xml.etree.ElementTree as ET


def get_mask(path):
        """
        Read values from RunInfo.xml needed to create a mask
        :param path: Full path to RunInfo.xml
        :return: Mask read from RunInfo.xml. A triplet of: Number (Read number), NumCycles, IsIndexedRead
        """

        # List of information required to be used as the mask in bcl2fastq
        # Groups of three elements: (Read) Number, NumCyles and IsIndexedRead will be stored consecutively
        mask = []
        file_name = path + '/RunInfo.xml'
        # get the tree of the XML file
        tree = ET.parse(file_name).getroot()
        # we are only intrested in the the Reads section
        root = tree.find('Run/Reads').getchildren()
        # loop over all child elements of Reads storing (in order): (read) Number, NumCycles, and IsIndexedRead
        for i in root:
            mask.append(i.get('Number'))
            mask.append(i.get('NumCycles'))
            mask.append(i.get('IsIndexedRead'))

        return mask

# Unit Test
if __name__ == '__main__':
    file_path = "/home/U008/lcebaman/scripts/data/RunInfo.xml"
    print(file_path)

    mask = get_mask(file_path)
    print(mask)
