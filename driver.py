#!/opt/anaconda/bin/python

# Runs BCL2fastq and then fastqc on a given run directory

import os
import logging
import argparse
from datetime import datetime

from utils import xmlparsing
from utils import sample_sheet_parser
from utils import create_bcl2fastq_PBS
from utils import create_fastqc_PBS
from utils import qsub_dependents

import config


if __name__ == '__main__':

    # parse the input directory /abs/path/to/INPUT_DATA/run_name
    parser = argparse.ArgumentParser()
    parser.add_argument('dirname')
    args = parser.parse_args()

    run_name = os.path.basename(args.dirname) 
    fastq_path = os.path.join(config.fastq, run_name) + '/'
    job_dir = os.path.join(config.jobs, run_name) + '/'
    
    if not os.path.exists(fastq_path):
        os.makedirs(fastq_path)
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    log_file = job_dir + 'driver_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log'
    logging.basicConfig(
        filename=log_file, format='%(asctime)s %(message)s', datefmt='[%d/%m/%Y-%H:%M:%S]', level=logging.INFO)
    logger = logging.getLogger('AnalysisDriver')
    logger.info('Reading bcl data from %s ', args.dirname)
    logger.info('Fastq path is %s ', fastq_path)

    # Read RunInfo.xml
    logger.info('Reading the reads info from ' + args.dirname)
    read_details = xmlparsing.get_mask(args.dirname)
    logger.info('Read details [ReadNumber, NumCycles, IsIndexedRead, ...] are %s ',read_details)
    
    # Create BCL2FASTQ PBS script
    logger.info('Create BCL2FASTQ pbs script')
    bcl2fastq_PBS_name = 'BCL2FASTQ_' + run_name + '.pbs'
    logger.info('bcl2fastq PBS File is ' + bcl2fastq_PBS_name)

    create_bcl2fastq_PBS.bcl2fastq_PBS(read_details, job_dir+bcl2fastq_PBS_name, args.dirname, fastq_path)

    # Create the fastqc PBS script
    logger.info('Create fastqc PBS script')
    fastqc_PBS_name = 'FASTQC_' + run_name + '.pbs'
    logger.info('fastqc PBS File is ' + fastqc_PBS_name)

    create_fastqc_PBS.fastqc_PBS(job_dir + fastqc_PBS_name, fastq_path)

    os.chdir(job_dir)
    
    # submit the BCL2FASTQ script to batch scheduler
    logger.info('Submitting: ' + job_dir + bcl2fastq_PBS_name)
    bcl2fastq_jobid = qsub_dependents.qsub([job_dir+bcl2fastq_PBS_name])
    logger.info('BCL2FASTQ jobId: ' + bcl2fastq_jobid)
    
    # submit the fastqc scrpipt to the batch scheduler
    logger.info('Submitting: %s', job_dir+fastqc_PBS_name)
    jobid = qsub_dependents.qsub_dependents([job_dir+fastqc_PBS_name], jobid=bcl2fastq_jobid)
    logger.info('FASTQC jobId: %s', jobid)
    logger.info('Done')

