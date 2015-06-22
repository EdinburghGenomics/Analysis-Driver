def get_mask(file_path):
    """
    Reads the values needed to create the mask.
    :param file_path: Full path to RunInfo.xml
    :return: A mask read from RunInfo.xml. The mask is made of trio values which indicate:
        1 - Number (the read number)
        2 - NumCycles
        3 - IsIndexedRead
    """

    import xml.etree.ElementTree as ET

    # List of information required to be used as the mask in bcl2fastq
    # Groups of three elements: (Read) Number, NumCyles and IsIndexedRead will be stored consecutively
    mask = []
    tree = ET.parse(file_path + '/RunInfo.xml').getroot()
    element = tree.find('Run/Reads')  # we are only interested in the the Reads section

    # Loop over all child elements of Reads storing (in order): (read) Number, NumCycles, and IsIndexedRead
    for i in element.getchildren():
        mask.append(i.get('Number'))
        mask.append(i.get('NumCycles'))
        mask.append(i.get('IsIndexedRead'))
     
    return mask


# Unit test
if __name__ == '__main__':
    file_path = '/home/U008/lcebaman/scripts/data/RunInfo.xml'
    print(file_path)
    mask = get_mask(file_path)
    print(mask)
