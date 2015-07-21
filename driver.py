import os
import logging
import argparse
from datetime import datetime
from time import sleep

import pbs_executor
from utils import xmlparsing
from utils import sample_sheet_parser
from utils import qsub_dependents
import util

import config as cfg

def main():

    # parse the input directory /abs/path/to/INPUT_DATA/run_name
    parser = argparse.ArgumentParser()
    parser.add_argument('dirname')
    args = parser.parse_args()

    config = cfg.Configuration()

    run_name = os.path.basename(args.dirname) 
    fastq_path = os.path.join(config.fastq, run_name)
    job_dir = os.path.join(config.jobs, run_name)
    
    if not os.path.exists(fastq_path):
        os.makedirs(fastq_path)
    if not os.path.exists(job_dir):
        os.makedirs(job_dir)

    log_file = os.path.join(job_dir, 'driver_' + datetime.now().strftime('%Y%m%d-%H%M%S') + '.log')
    logging.basicConfig(
        # stream=sys.stdout,
        stream=open(log_file, 'w'),
        format='%(asctime)s %(message)s',
        datefmt='[%d-%m-%Y %H:%M:%S]',
        level=logging.INFO
    )
    logger = logging.getLogger('AnalysisDriver')
    logger.info('Reading bcl data from ' + args.dirname)
    logger.info('Fastq path is ' + fastq_path)

    # Read RunInfo.xml
    logger.info('Reading the reads info from ' + args.dirname)

    sample_sheet = sample_sheet_parser.SampleSheet(args.dirname)
    run_info = xmlparsing.RunInfo(args.dirname)
    if run_info.barcode_len:
        if not sample_sheet.check_barcodes() == run_info.barcode_len:
            logger.warn('Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' % (
                sample_sheet.check_barcodes(), run_info.barcode_len
                )
            )
    else:
        logger.warn('No barcode in RunInfo.xml')

    mask = run_info.mask.tostring(sample_sheet.check_barcodes())
    # mask = 'y150n,i6,y150n'
    logger.info('bcl2fastq mask: ' + mask)
    
    # Create BCL2FASTQ PBS script
    logger.info('Create bcl2fastq PBS script')
    bcl2fastq_pbs_name = os.path.join(job_dir, 'bcl2fastq_' + run_name + '.pbs')
    logger.info('bcl2fastq PBS file is ' + bcl2fastq_pbs_name)
    bcl2fastq_writer = pbs_executor.BCL2FastqPBSWriter(
        bcl2fastq_pbs_name, 'bcl2fastq', os.path.join(job_dir, 'bcl2fastq_pbs.log')
    )
    bcl2fastq_writer.write(mask, args.dirname, fastq_path)

    # Create the fastqc PBS script
    logger.info('Creating fastqc PBS script')
    fastqc_pbs_name = os.path.join(job_dir, 'fastqc_' + run_name + '.pbs')
    logger.info('Fastqc PBS file is ' + fastqc_pbs_name)

    fastqc_writer = pbs_executor.FastqcPBSWriter(
        fastqc_pbs_name, 'fastqc', os.path.join(job_dir, 'fastqc_pbs.log')
    )
    fastqc_writer.write(fastq_path, job_dir)

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

        fastqc_complete = os.path.join(job_dir, '.fastqc_complete')
        logger.info('Waiting for creation of ' + fastqc_complete)
        while not os.path.exists(fastqc_complete):
            sleep(15)

        logger.info('Fastqc complete, executing alignment')
    
    csv_writer = pbs_executor.BCBioCSVWriter(fastq_path, job_dir, sample_sheet)
    csv_writer.write()
    fastqs = []

    for sample_project in sample_sheet.sample_projects:
        fastqs = fastqs + util.find_fastqs(os.path.join(fastq_path, sample_project))

    os.chdir(job_dir)

    bcbio_pbs = os.path.join(job_dir, 'bcbio.pbs')
    bcbio_writer = pbs_executor.BCBioPBSWriter(bcbio_pbs, 'bcbio', 'bcbio.log')
    bcbio_writer.setup_bcbio_run(
        config.bcbio,
        os.path.join(os.path.dirname(__file__), 'etc', 'bcbio_alignment.yaml'),
        os.path.join(job_dir, 'bcbio'),
        os.path.join(job_dir, 'samples.csv'),
        fastqs
    )

    bcbio_writer.write(
        config.bcbio,
        os.path.join(job_dir, 'samples', 'config', 'samples.yaml'),
        os.path.join(job_dir, 'samples', 'work')
    )
    bcbio_jobid = str(
        # qsub_dependents.qsub_dependents([bcbio_pbs], jobid=fastqc_jobid)
        qsub_dependents.qsub([bcbio_pbs])
    ).lstrip('b\'').rstrip('\'')
    logger.info('BCBio jobId: ' + bcbio_jobid)

    logger.info('Done')

if __name__ == '__main__':
    main()

