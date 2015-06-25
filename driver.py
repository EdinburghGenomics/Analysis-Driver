#!/opt/anaconda/bin/python

# Runs BCL2fastq and then fastqc on a given run directory

import sys, getopt
import os,logging
import argparse
sys.path.append('imports')

from xmlparsing import getMask
from csvparsing import readSampleSheet,getSampleProject
from createBCL2FASTQ_PBS import generateMask,bcl2fastq_PBS
from createFASTQC_PBS import fastqc_PBS
from subprocess import call
from qsub_dependents import qsub,qsub_dependents
from datetime import datetime

if __name__ == "__main__":

    # parse the input directory /abs/path/to/INPUT_DATA/runname
    parser = argparse.ArgumentParser()
    parser.add_argument('dirname')
    cmdargs = parser.parse_args()

    # We keep all paths with trailing '/' remove and add to ensure it is there
    inputPath = cmdargs.dirname.rstrip('/')+'/'
    runName   = os.path.basename(inputPath.rstrip('/')) 
    fastqPath = os.path.normpath(inputPath+'../../fastq/'+runName)+'/'
    jobsPath  = os.path.normpath(inputPath+'../../jobs/'+runName)+'/'

    # create project directory
    if not os.path.exists(fastqPath): os.makedirs(fastqPath)
    if not os.path.exists(jobsPath): os.makedirs(jobsPath)


    #open logging file and configure it with date and time stamp
    logfileName = jobsPath+'driver_'+datetime.now().strftime("%Y%m%d-%H%M%S")+'.log'
    logging.basicConfig(filename=logfileName,format='%(asctime)s %(message)s',
                        datefmt='[%d/%m/%Y-%H:%M:%S]',level=logging.INFO)
    logging.info('Reading bcl data from %s ',inputPath)
    logging.info('Fastq path is %s ',fastqPath)


    # Read RunInfo.xml
    logging.info('Reading the reads info from %s ',inputPath)
    readDetails = getMask(inputPath)
    logging.info('Read details [ReadNumber, NumCycles, IsIndexedRead, ...] are %s ',readDetails)
    
    # Create BCL2FASTQ PBS script
    logging.info('Create BCL2FASTQ pbs script')
    bcl2fastqPbsName = 'BCL2FASTQ_' + runName + '.pbs'
    logging.info('bcl2fastq PBS File is %s ',bcl2fastqPbsName)

    bcl2fastq_PBS(readDetails, jobsPath+bcl2fastqPbsName, inputPath, fastqPath)

    # Create the fastqc PBS script
    logging.info('Create fastqc PBS script')
    fastqcPbsName = 'FASTQC_' +runName + '.pbs'
    logging.info('fastqc PBS File is %s ' , fastqcPbsName)

    fastqc_PBS(jobsPath+fastqcPbsName, fastqPath)

    # get into pbs directory
    os.chdir(jobsPath)
    
    # submit the BCL2FASTQ script to batch scheduler
    logging.info('Submitting: %s',jobsPath+bcl2fastqPbsName)
    BCL2FASTQ_jobid = qsub([jobsPath+bcl2fastqPbsName])
    logging.info('BCL2FASTQ jobId: %s',BCL2FASTQ_jobid)
    
    # submit the fastqc scrpipt to the batch scheduler
    logging.info('Submitting: %s',jobsPath+fastqcPbsName)
    jobid = qsub_dependents([jobsPath+fastqcPbsName], jobid=BCL2FASTQ_jobid)
    logging.info('FASTQC jobId: %s',jobid)
