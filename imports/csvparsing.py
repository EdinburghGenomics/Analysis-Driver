
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
     d = {}
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
     # get the position in the CSV file (required for file naming)
     firstTime={}
     for row in dlist:
          if row[1] not in firstTime:
               firstTime[row[1]] = row[0]

# Store keys (sampleID) and values (Lane,SampleName, firstTimePos, projectName)

# {'10015AT0001L01': [['1', '50293', '1', '10015AT']],
#  '10015AT0002L01': [['2', '128864', '2', '10015AT']],
#  '10015AT0003L01': [['3', '172937', '3', '10015AT']],
#  '10015AT0004L01': [['4', 'Na12878', '4', '10015AT'],
#   ['6', 'Na12878', '4', '10015AT'],
#   ['7', 'Na12878', '4', '10015AT']],
#  '10015AT0005L01': [['5', 'PhiX', '5', '10015AT'],
#   ['8', 'PhiX', '5', '10015AT']]}


     for row in dlist:
          d.setdefault(row[1],[]).append([row[0],row[2],firstTime[row[1]],row[7]] )
     
     return len(dlist),d

# gives the sample Project
def getSampleProject(d):
     #the first one for instance
     return d.values()[0][0][3]
                                             
