#!/opt/anaconda/bin/python

#This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys, getopt
import os,csv
import subprocess
from collections import defaultdict
import itertools
from itertools import groupby
from operator import itemgetter

sys.path.append('imports')

from args import parsArgs

def getId(id):
    return id



if __name__ == "__main__":
   
    # parse the input directory
    inputDirectory = parsArgs(sys.argv[1:])

    # print inputDirectory
    d = defaultdict(list)
    reader = csv.reader(open(sys.argv[1], "rb"))
    
    # read lines until [Data] marker
    while not next(reader)[0].startswith('[Data]'):
        None
        next(reader) # ignore names

    #get rid of labels
    next(reader)
    for k, g in groupby(zip(reader, itertools.count()), key=lambda x: x[0][1]):
        map(lambda x: d[k].append((x[0][0], x[1], x[0][2])), g)

    FASTQ="NULL"  
    BASE_PATH= inputDirectory+"/Data/Intensities/BaseCalls"
    counter = 0   
    
    for k, v in sorted(d.items(), key=itemgetter(0)):
        #print 'id: %s, lanes: %s' % (k, v)
        FASTQ = BASE_PATH+"/"+str(v[0][2])+"/"+str(k)
        print len(v)
        for i in xrange(0,len(v)):
            
            if v[0][2] != '':
                FASTQ = BASE_PATH+"/"+str(v[0][2])+"/"+str(k)
                FASTQ1 =FASTQ+"/"+str(v[0][2])+"_S"+str(i+1)+"_L00"+str(v[i][0])+"_R"+str(i+1)+"_001"
                print FASTQ1
    #    else:
     #       FASTQ = BASE_PATH+"/"+k
           
      

        


