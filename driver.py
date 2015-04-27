#!/opt/anaconda/bin/python

#This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys, getopt
import os,logging
sys.path.append('imports')

from xmlparsing import getMask
from csvparsing import readSampleSheet
from createBCL2FASTQ_PBS import generateMask,bcl2fastq_PBS
from createBCBIO_PBS import bcbio_loop
from subprocess import call
from qsub_dependents import qsub,qsub_dependents
from args import parsArgs
from datetime import datetime
 

if __name__ == "__main__":

    #open logging file and configure it with date and time stamp
    logfileName = 'log_'+datetime.now().strftime("%Y%m%d-%H%M%S")+'.log'
    logging.basicConfig(filename=logfileName,format='%(asctime)s %(message)s', datefmt='[%d/%m/%Y-%H:%M:%S]')

    # parse the input directory
    inputDirectory = parsArgs(sys.argv[1:])
    
    logging.warning('Reading bcl data from %s ',inputDirectory)
    
    # create pbs temporary directory
    if not os.path.exists("pbs"):
        os.makedirs("pbs")

    
    # Read RunInfo.xml
    logging.warning('Reading the mask from %s ',inputDirectory)
    mask = getMask(inputDirectory)
    
    # Read SampleSheet.csv
    logging.warning('Reading SampleSheet from %s ',inputDirectory)
    sheetDict = readSampleSheet(inputDirectory)
    
    # Create BCL2FASTQ PBS script
    logging.warning('Create BCL2FASTQ pbs script')
    bcl2fastq_PBS(mask, inputDirectory)

    # create BCBIO PBS scripts
    logging.warning('Creating BCBIO PBS scripts')
    bcbio_loop(sheetDict, inputDirectory)
    
    # get into pbs directory
    os.chdir("pbs")
    
    # submit bcl2fastq 
    # create a list with the name of the PBS script
    argsString='BCL_'+ str(inputDirectory)+'.pbs'
    args=[argString]
    
    # submit the BCL2FASTQ script to bach scheduler
    logging.warning('Submit BCL2FASTQ_PBS')
    #BCL2FASTQ_jobid = qsub_dependents(args)
    
    # submit the BCBIO scripts once BCL2FASTQ has finished
    logging.warning('Submit BCBIO_PBS')

#     list_jobIds=[]
#     # submit set of BCBIO jobs. A job per sampleId included in the SampleSheet
#     for i in range(0, len(d)):
#        scriptName = "runBCBIO_" +`i`+".pbs"
#        args=[scriptName]
#        # store the jobIds in a list. They will not get executed until BCL2FASTQ has finished
#        list_jobIds.append( qsub_dependents(args,BCL2FASTQ_jobid))

    logging.warning('Submit BCBIO_PBS')
    
