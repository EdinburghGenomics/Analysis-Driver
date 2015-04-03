
# reads the sampleIds and sampleNames 
# Input: Full path to SampleSheet.csv
# Outout: It returns to lists
#         1 - list of sampleIds
#         2-  list of sampleNames


def readSampleSheet(filename):
     
     import os,sys
     import codecs,csv

     csvlist=[]
     sampleId=[]
     sampleName=[]
     sampleProject=[]


     # read the whole CSV file
     try:
          with codecs.open(filename, "rb") as f_obj:
               reader = csv.reader(f_obj)
               for row in reader:
                    csvlist.append(row)
                    
     except IOError:
          print "No such file : ", filename
          print "Please, introduce a valid path to SampleSheet.csv file"
          sys.exit()
          
      #find [Data] tag
     dataTagPos = 0
     for i in xrange(10,len(csvlist)):
          if csvlist[i][0] == '[Data]':
               dataTagPos = i
                 
       #First row of reald data
     dataPos=dataTagPos + 2

       # store values in lists
     for i in xrange(dataPos,len(csvlist)):
          sampleId.append(csvlist[i][1])
          sampleName.append(csvlist[i][2])
          sampleProject.append(csvlist[i][7])

     return sampleId,sampleName


# Unit Test
filename = "/home/U008/lcebaman/scripts/data/SampleSheet.csv"
print filename
sampleId,sampleName = readSampleSheet(filename)

