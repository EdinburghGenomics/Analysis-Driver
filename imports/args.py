def parsArgs(argv):
   import sys, getopt
   import os

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

   print 'Input directory is', inputfile.rstrip('/')
   return inputfile.rstrip('/')
