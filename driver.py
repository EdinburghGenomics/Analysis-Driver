#!/opt/anaconda/bin/python

#This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys, getopt
import os
from xmlparsing import getMask
from csvparsing import readSampleSheet
from createBCL2FASTQ_PBS import generateMask,bcl2fastq_PBS
def parsArgs(argv):
   inputfile = ''
 
   try:
      opts, args = getopt.getopt(argv,"hi:",["ifile="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> '
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> '
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
         assert os.path.exists(inputfile), "Non existing directory, "+ str(inputfile)
      else:
           print 'test.py -i <inputfile> '
           sys.exit()

   print 'Input directory is', inputfile
   return inputfile

if __name__ == "__main__":
   # parse the input directory
   inputDirectory = parsArgs(sys.argv[1:])
   #print inputDirectory
   
   # Read RunInfo.xml
   mask = getMask(inputDirectory)
   # print mask

   # Read SampleSheet.csv
   sampleId,sampleName = readSampleSheet(inputDirectory)
   # print sampleId,sampleName

   # Create BCL2FASTQ PBS script
   bcl2fastq_PBS(mask)
