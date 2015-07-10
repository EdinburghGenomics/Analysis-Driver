#!/opt/anaconda/bin/python

# Runs BCL2fastq and then fastqc on a given run directory

import os
import logging
import argparse
from datetime import datetime
import subprocess

from utils import xmlparsing
from utils import sample_sheet_parser
import pbs_executor

from utils import qsub_dependents
import util
import config


def main():

    sample_id = '10015AT'

    # parse the input directory /abs/path/to/INPUT_DATA/run_name
    parser = argparse.ArgumentParser()
    parser.add_argument('dirname')
    args = parser.parse_args()

    run_name = os.path.basename(args.dirname) 
    fastq_path = os.path.join(config.fastq, run_name)
    job_dir = os.path.join(config.jobs, run_name)
    
    if not os.path.exists(fastq_path):
        os.makedirs(fastq_path)
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    log_file = os.path.join(job_dir, 'driver_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log')
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
    mask = xmlparsing.get_mask(args.dirname)
    logger.info('Read details [ReadNumber, NumCycles, IsIndexedRead, ...] are ' + str(mask))
    
    # Create BCL2FASTQ PBS script
    logger.info('Create BCL2FASTQ pbs script')
    bcl2fastq_pbs_name = os.path.join(job_dir, 'bcl2fastq_' + run_name + '.pbs')
    logger.info('bcl2fastq PBS File is ' + bcl2fastq_pbs_name)
    bcl2fastq_writer = pbs_executor.BCL2FastqPBSWriter(
        bcl2fastq_pbs_name, 'bcl2fastq', os.path.join(job_dir, 'bcl2fastq_pbs.log')
    )
    bcl2fastq_writer.write(mask, args.dirname, fastq_path)

    # Create the fastqc PBS script
    logger.info('Create fastqc PBS script')
    fastqc_pbs_name = os.path.join(job_dir, 'fastqc_' + run_name + '.pbs')
    logger.info('fastqc PBS File is ' + fastqc_pbs_name)

    fastqc_writer = pbs_executor.FastqcPBSWriter(
        fastqc_pbs_name, 'fastqc', os.path.join(job_dir, 'fastqc_pbs.log')
    )
    fastqc_writer.write(fastq_path)

    sample_projects = [os.path.abspath(proj) for proj in os.listdir(os.path.join(fastq_path, sample_id))]

    bcbio_jobs = sample_sheet_parser.read_sample_sheet(args.dirname)
    # {
    #     sample_id: (lane_7, sample_name, position),
    #     sample_id: (lane_7, sample_name, position)
    # }

    # Write the bcbio PBS script
    for sample_project in sample_projects:
        bcbio_run_dir = os.path.join(job_dir, 'bcbio', sample_project)
        if not os.path.exists(bcbio_run_dir):
            os.mkdir(bcbio_run_dir)

        bcbio_pbs = os.path.join(job_dir, 'bcbio_' + os.path.basename(sample_project) + '.pbs')
        bcbio_writer = pbs_executor.BCBioPBSWriter(
            bcbio_pbs, 'bcbio_alignment', os.path.join(bcbio_run_dir, 'log.txt')
        )
        fastqs = bcbio_writer.get_fastqs(fastq_path, run_name, sample_id, sample_project)
        bcbio_writer.setup_bcbio_run(
            config.bcbio,
            os.path.join(
                os.path.dirname(__file__), 'etc', 'bcbio_alignment.yaml'
            ),
            sample_project,
            fastqs
        )
        bcbio_writer.write(
            config.bcbio,
            os.path.join(job_dir, 'bcbio', sample_project + '.yaml'),
            os.path.join(job_dir, 'bcbio', sample_project, 'work')
        )

    bcbio_pbs_scripts = [
        os.path.abspath(f) for f in os.listdir(job_dir) if 'bcbio' in f and f.endswith('.pbs')
    ]
    # submit the BCL2FASTQ script to batch scheduler
    if config.job_execution == 'pbs':
        logger.info('Submitting: ' + os.path.join(job_dir, bcl2fastq_pbs_name))
        bcl2fastq_jobid = str(
            qsub_dependents.qsub([os.path.join(job_dir, bcl2fastq_pbs_name)])
        ).lstrip('b\'').rstrip('\'')
        logger.info('BCL2FASTQ jobId: ' + bcl2fastq_jobid)
    
        # submit the fastqc script to the batch scheduler
        logger.info('Submitting: ' + os.path.join(job_dir, fastqc_pbs_name))
        fastqc_jobid = str(
            qsub_dependents.qsub_dependents([os.path.join(job_dir, fastqc_pbs_name)], jobid=bcl2fastq_jobid)
        ).lstrip('b\'').rstrip('\'')
        logger.info('FASTQC jobId: ' + fastqc_jobid)

        # Submit the bcbio PBS scripts
        for script in bcbio_pbs_scripts:
            logger.info('Submitting bcbio job: ' + script)
            bcbio_jobid = str(
                qsub_dependents.qsub_dependents([script], jobid=fastqc_jobid)
            ).lstrip('b\'').rstrip('\'')
            logger.info('BCBIO jobId: ' + bcbio_jobid)

    elif config.job_execution == 'local':
        logger.info('Executing locally: ' + os.path.join(job_dir, bcl2fastq_pbs_name))
        util.localexecute('sh', os.path.join(job_dir, bcl2fastq_pbs_name))

        logger.info('Executing locally: ' + os.path.join(job_dir, fastqc_pbs_name))
        util.localexecute('sh', os.path.join(job_dir, bcl2fastq_pbs_name))

        for script in bcbio_pbs_scripts:
            logger.info('Executing locally: ' + script)
            util.localexecute('sh', script)

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

