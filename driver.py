#!/opt/anaconda/bin/python

#This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys, getopt
import os
sys.path.append('imports')

from xmlparsing import getMask
from csvparsing import readSampleSheet
from createBCL2FASTQ_PBS import generateMask,bcl2fastq_PBS
from createBCBIO_PBS import bcbio_loop
from subprocess import call
from qsub_dep import qsub,qsub_dependents
from args import parsArgs


if __name__ == "__main__":
   
    # parse the input directory
    inputDirectory = parsArgs(sys.argv[1:])
    # print inputDirectory
    
    # create pbs temporary directory
    if not os.path.exists("pbs"):
        os.makedirs("pbs")
    
    # Read RunInfo.xml
    mask = getMask(inputDirectory)
    # print mask

    # Read SampleSheet.csv
    sampleId,sampleName = readSampleSheet(inputDirectory)
    # print sampleId,sampleName

    # Create BCL2FASTQ PBS script
    bcl2fastq_PBS(mask)

    # Total number of runs, one per sampleId
    
    # create BCBIO PBS scripts
    bcbio_loop(sampleId, sampleName, inputDirectory)
    
    # get into pbs directory
    os.chdir("pbs")
    # submit bcl2fastq 
    # create a list with the name of the PBS script
    args=['runBCL2FASTQ.pbs']

    # submit the BCL2FASTQ script to bach scheduler
    BCL2FASTQ_jobid = qsub_dependents(args)

    list_jobIds=[]
    # submit set of BCBIO jobs
    for i in range(0, n):
       scriptName = "runBCBIO" +`i`+".pbs"
       args=[scriptName]
       # store the jobIds in a list
       list_jobIds.append( qsub_dependents(args,BCL2FASTQ_jobid))

    
    # TODO: once jobs are done, we need to move the data output back to RDF
    #       Need to discuss how we can do it. Either from the same job or submitting new jobs
    #       Also the resources allocated for copying back
