
# reads the sampleIds and sampleNames 
# Input: Full path to SampleSheet.csv
# Outout: It returns to lists
#         1 - list of sampleIds
#         2-  list of sampleNames

def readSampleSheet(file_path):
     import subprocess
     from collections import defaultdict
     import itertools
     from itertools import groupby
     from operator import itemgetter
     import csv

     # create a dictionary
     d = defaultdict(list)
     # temp list
     dlist= []
     # hardcoded name 
     filename = file_path + "/SampleSheet.csv"
     # open the file with csv
     reader = csv.reader(open(filename, "rb"))
    
     # read lines until [Data] marker
     while not next(reader)[0].startswith('[Data]'):
          None
     next(reader) # ignore names

     # store the rest of the csv file in a list
     for row in reader:
          # store only the non empty lines
          if any(row):
               dlist.append(row)
     
     # print len(dlist)
     # loop = len(dlist) -1
     # print "loop= ",loop
     # #make sure there are no empyt lines
     # for i in xrange(0,loop):
     #      print dlist[i]
     #      if not dlist[i]:
     #           del(dlist[i])
     
     print len(dlist)
     # stores elements as key(sampleId) -> values (Lane, Position, sampleName)
     for k, g in groupby(zip(dlist, itertools.count()), key=lambda x: x[0][1]):
          map(lambda x: d[k].append((x[0][0], x[1], x[0][2])), g)

     return d

