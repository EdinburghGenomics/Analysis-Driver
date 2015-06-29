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
from utils import create_bcbio_PBS
from utils import qsub_dependents

send_qsubs = True


if __name__ == '__main__':

    # parse the input directory /abs/path/to/INPUT_DATA/run_name
    parser = argparse.ArgumentParser()
    parser.add_argument('dirname')
    args = parser.parse_args()

    run_name = os.path.basename(args.dirname) 
    fastq_path = os.path.normpath(os.path.join(args.dirname, '..', '..', 'fastq', run_name)) + '/'
    job_dir = os.path.normpath(os.path.join(args.dirname, '..', '..', 'jobs', run_name)) + '/'
    
    if not os.path.exists(fastq_path):
        os.makedirs(fastq_path)
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    log_file = job_dir + 'driver_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log'
    logging.basicConfig(
        filename=log_file,
        format='%(asctime)s %(message)s',
        datefmt='[%d/%m/%Y-%H:%M:%S]',
        level=logging.INFO
    )
    logger = logging.getLogger('AnalysisDriver')
    logger.info('Reading bcl data from ' + args.dirname)
    logger.info('Fastq path is ' + fastq_path)

    # Read RunInfo.xml
    logger.info('Reading the reads info from ' + args.dirname)
    read_details = xmlparsing.get_mask(args.dirname)
    logger.info('Read details [ReadNumber, NumCycles, IsIndexedRead, ...] are ' + str(read_details))
    
    # Create BCL2FASTQ PBS script
    logger.info('Create BCL2FASTQ pbs script')
    bcl2fastq_PBS_name = 'BCL2FASTQ_' + run_name + '.pbs'
    logger.info('bcl2fastq PBS File is ' + bcl2fastq_PBS_name)
    create_bcl2fastq_PBS.bcl2fastq_PBS(
        read_details,
        job_dir + bcl2fastq_PBS_name,
        args.dirname,
        fastq_path
    )

    # Create the fastqc PBS script
    logger.info('Create fastqc PBS script')
    fastqc_PBS_name = 'FASTQC_' + run_name + '.pbs'
    logger.info('fastqc PBS File is ' + fastqc_PBS_name)
    create_fastqc_PBS.fastqc_PBS(job_dir + fastqc_PBS_name, fastq_path)

    # Write the bcbio PBS script
    logger.info('Writing bcbio script')
    bcbio_PBS_name = 'BCBIO_' + run_name
    logger.info('bcbio file: ' + bcbio_PBS_name)

    bcbio_jobs = sample_sheet_parser.read_sample_sheet(args.dirname)
    # {
    #     sample_id: (lane_7, sample_name, position),
    #     sample_id: (lane_7, sample_name, position)
    # }

    # pbs_name, bcbio_run_folder, run_id, lanes, sample_name='Unassigned_S0'
    bcbio_lanes = [x[0] for x in bcbio_jobs.values()[0]]
    for lane in bcbio_lanes:
        create_bcbio_PBS.bcbio_PBS(
            job_dir + bcbio_PBS_name,
            job_dir + 'bcbio/',
            run_name,
            lane  # ,
            # sample_sheet_parser.get_sample_project(bcbio_jobs)
        )

    os.chdir(job_dir)
    
    # submit the BCL2FASTQ script to batch scheduler
    logger.info('Submitting: ' + job_dir + bcl2fastq_PBS_name)
    if send_qsubs:
        bcl2fastq_jobid = str(
            qsub_dependents.qsub([job_dir + bcl2fastq_PBS_name])
        ).lstrip('b\'').rstrip('\'')
        logger.info('BCL2FASTQ jobId: ' + bcl2fastq_jobid)
    
    # submit the fastqc scrpipt to the batch scheduler
    logger.info('Submitting: ' + job_dir + fastqc_PBS_name)
    if send_qsubs:
        fastqc_jobid = str(
            qsub_dependents.qsub_dependents([job_dir + fastqc_PBS_name], jobid=bcl2fastq_jobid)
        ).lstrip('b\'').rstrip('\'')
        logger.info('FASTQC jobId: ' + fastqc_jobid)

    # Submit the bcbio PBS script
    logger.info('Submitting bcbio script')
    if send_qsubs:
        for lane in bcbio_lanes:
            bcbio_jobid = str(
                qsub_dependents.qsub_dependents(
                    [job_dir + bcbio_PBS_name + '_L00' + lane + '.pbs'],
                    jobid=fastqc_jobid
                )
            ).lstrip('b\'').rstrip('\'')
            logger.info('BCBIO jobId: ' + bcbio_jobid)

    logger.info('Done')

