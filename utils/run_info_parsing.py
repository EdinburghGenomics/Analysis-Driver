import xml.etree.ElementTree as ET


def get_mask(path):
    """
    Read the values needed to create a mask
    :param path: Full path to RunInfo.xml
    :return: A list containing: the read number, the number of cycles and whether the read is indexed
    """

    # List of information required to be used as the mask in bcl2fastq
    # Groups of three elements: (Read) Number, NumCyles and IsIndexedRead will be stored consecutively
    mask = []
    filename = path + '/RunInfo.xml'
    # get the tree of the XML file
    root = ET.parse(filename).getroot()
    # we are only interested in the the Reads section
    reads = root.find('Run/Reads').getchildren()

    # loop over all child elements of Reads storing (in order): (read) Number, NumCycles, and IsIndexedRead
    for i in reads:
        mask.append(i.get('Number'))
        mask.append(i.get('NumCycles'))
        mask.append(i.get('IsIndexedRead'))

    return mask


if __name__ == '__main__':
    file_path = '/home/U008/lcebaman/scripts/data/RunInfo.xml'
    print(file_path)
    mask = get_mask(file_path)
    print(mask)
