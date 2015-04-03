

# reads the values needed to create the mask.
# Input: Full path to RunInfo.xml
# Outout: Mask read from RunInfo.xml. The mask is made of pair values which indicate:
#         1 - NumCycles
#         2-  IsIndexedRead

def getMask(filename):

     import xml.etree.ElementTree as ET
     import os,sys

     # List of information required to be used as the mask in bcl2fastq
     # Groups of two elements: NumCyles and IsIndexedRead will be stored consecutively
     mask=[]
     
     # get the tree of the XML file
     tree = ET.parse(file_path).getroot()
     # we are only intrested in the the Reads section
     object = tree.find('Run/Reads')
     #loop over all child elements of Reads storing(in order) NUmCycles and IsIndexedRead info
     for i in  object.getchildren():
          mask.append(i.get('NumCycles'))
          mask.append(i.get('IsIndexedRead'))
     
     return mask

# Unit Test
file_path = "/home/U008/lcebaman/scripts/data/RunInfo.xml"
print file_path

mask = getMask(file_path)

print mask
