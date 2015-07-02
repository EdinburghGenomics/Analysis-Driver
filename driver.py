#!/opt/anaconda/bin/python

# Runs BCL2fastq and then fastqc on a given run directory

import os
import logging
import argparse
from datetime import datetime
import subprocess

from utils import xmlparsing
from utils import sample_sheet_parser
from utils import create_bcl2fastq_PBS
from utils import create_fastqc_PBS
from utils import create_bcbio_PBS
from utils import qsub_dependents

import config

job_execution = 'pbs'  # 'pbs' 'local' None


def main():

    # parse the input directory /abs/path/to/INPUT_DATA/run_name
    parser = argparse.ArgumentParser()
    parser.add_argument('dirname')
    args = parser.parse_args()

    run_name = os.path.basename(args.dirname) 
    fastq_path = os.path.join(config.fastq, run_name) + '/'
    job_dir = os.path.join(config.jobs, run_name) + '/'
    
    if not os.path.exists(fastq_path):
        os.makedirs(fastq_path)
    if not os.path.exists(job_path):
        os.makedirs(job_path)

    log_file = os.path.join(job_path, 'driver_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log')
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
        os.path.join(job_path, bcl2fastq_PBS_name),
        args.dirname,
        fastq_path
    )

    # Create the fastqc PBS script
    logger.info('Create fastqc PBS script')
    fastqc_PBS_name = 'FASTQC_' + run_name + '.pbs'
    logger.info('fastqc PBS File is ' + fastqc_PBS_name)
    create_fastqc_PBS.fastqc_PBS(os.path.join(job_path, fastqc_PBS_name), fastq_path)

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
    print(bcbio_jobs.values())
    bcbio_lanes = [x[0] for x in list(bcbio_jobs.values())[0]]
    for lane in bcbio_lanes:
        create_bcbio_PBS.bcbio_PBS(
            os.path.join(job_path, bcbio_PBS_name),
            os.path.join(job_path, 'bcbio'),
            run_name,
            lane,
            fastq_path  # ,
            # sample_sheet_parser.get_sample_project(bcbio_jobs)
        )

    os.chdir(job_path)
    
    # submit the BCL2FASTQ script to batch scheduler
    if job_execution == 'pbs':
        logger.info('Submitting: ' + os.path.join(job_path, bcl2fastq_PBS_name))
        bcl2fastq_jobid = str(
            qsub_dependents.qsub([os.path.join(job_path, bcl2fastq_PBS_name)])
        ).lstrip('b\'').rstrip('\'')
        logger.info('BCL2FASTQ jobId: ' + bcl2fastq_jobid)
    
        # submit the fastqc script to the batch scheduler
        logger.info('Submitting: ' + os.path.join(job_path, fastqc_PBS_name))
        fastqc_jobid = str(
            qsub_dependents.qsub_dependents([os.path.join(job_path, fastqc_PBS_name)], jobid=bcl2fastq_jobid)
        ).lstrip('b\'').rstrip('\'')
        logger.info('FASTQC jobId: ' + fastqc_jobid)

        # Submit the bcbio PBS script
        for lane in bcbio_lanes:
            logger.info('Submitting bcbio job for lane ' + lane)
            bcbio_jobid = str(
                qsub_dependents.qsub_dependents(
                    [os.path.join(job_path, bcbio_PBS_name + '_L00' + lane + '.pbs')],
                    jobid=fastqc_jobid
                )
            ).lstrip('b\'').rstrip('\'')
            logger.info('BCBIO jobId: ' + bcbio_jobid)

    elif job_execution == 'local':
        logger.info('Executing locally: ' + os.path.join(job_path, bcl2fastq_PBS_name))
        execute_bash(os.path.join(job_path, bcl2fastq_PBS_name))

        logger.info('Executing locally: ' + os.path.join(job_path, fastqc_PBS_name))
        execute_bash(os.path.join(job_path, bcl2fastq_PBS_name))

        for lane in bcbio_lanes:
            script = os.path.join(job_path, bcbio_PBS_name + '_L00' + lane + '.pbs')
            logger.info('Executing locally: ' + script)
            execute_bash(script)

    else:
        logger.info('No job_execution set. Scripts written but not executed.')

    logger.info('Done')


def execute_bash(script):
    proc = subprocess.Popen(['sh', script], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    print(out)
    print(err)


if __name__ == '__main__':
    main()

