#!/opt/anaconda/bin/python

#This program expects an input directory. From there, we will run BCL2FASTQ and BCBIO
import sys, getopt
import os,logging
sys.path.append('imports')

from xmlparsing import getMask
from csvparsing import readSampleSheet,getSampleProject
from createBCL2FASTQ_PBS import generateMask,bcl2fastq_PBS
from createBCBIO_PBS import bcbio_loop
from subprocess import call
from qsub_dependents import qsub,qsub_dependents
from args import parsArgs
from datetime import datetime
from makeProject import makeProject, getDirName

if __name__ == "__main__":

    #open logging file and configure it with date and time stamp
    logfileName = 'log_'+datetime.now().strftime("%Y%m%d-%H%M%S")+'.log'
    logging.basicConfig(filename=logfileName,format='%(asctime)s %(message)s', datefmt='[%d/%m/%Y-%H:%M:%S]')

    # parse the input directory
    inputPath = parsArgs(sys.argv[1:])
    
    logging.warning('Reading bcl data from %s ',inputPath)

    # create project directory
    projectName = getDirName(inputPath)  
    workDir = projectName+"_work"
    makeProject(workDir)

    # Read RunInfo.xml
    logging.warning('Reading the mask from %s ',inputPath)
    mask = getMask(inputPath)
    
    # Read SampleSheet.csv
    logging.warning('Reading SampleSheet from %s ',inputPath)
    numLanes,sheetDict = readSampleSheet(inputPath)
    print numLanes
    # get sampleProject
    sampleProject = getSampleProject(sheetDict)

    # Create BCL2FASTQ PBS script
    logging.warning('Create BCL2FASTQ pbs script')
    pbsName = 'BCL_'+ projectName +'.pbs'
    print pbsName

    bcl2fastq_PBS(mask, pbsName, workDir, inputPath)

    # create BCBIO PBS scripts
    logging.warning('Creating BCBIO PBS scripts')
    bcbio_loop(sheetDict, inputPath, workDir, sampleProject)
    
     # get into pbs directory
    os.chdir(workDir+"/pbs")
    
    # submit bcl2fastq 
    # create a list with the name of the PBS script

    args=[pbsName]
    
    # submit the BCL2FASTQ script to batch scheduler
    logging.warning('Submit BCL2FASTQ_PBS')
    BCL2FASTQ_jobid = qsub_dependents(args)
    
    # submit the BCBIO scripts once BCL2FASTQ has finished
    logging.warning('Submit BCBIO_PBS')

    list_jobIds=[]
     # submit set of BCBIO jobs. A job per sampleId included in the SampleSheet
    for i in range(1, numLanes):
        scriptName = "runBCBIO_" +`i`+".pbs"
        args=[scriptName]
         # store the jobIds in a list. They will not get executed until BCL2FASTQ has finished
        list_jobIds.append( qsub_dependents(args,BCL2FASTQ_jobid))
         
    logging.warning('Submit BCBIO_PBS')
        
