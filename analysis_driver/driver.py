import os
import sys
import argparse
import logging
import logging.config
from time import sleep

from analysis_driver import writer, util
from analysis_driver.legacy import qsub_dependents
import analysis_driver.reader
from analysis_driver.config import default as config  # imports the default config singleton

logging_helper = util.NamedAppLogger('driver')


def main():
    """
    :return: Exit status
    :rtype: int
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input_run_folder', type=str, help='An absolute path to an input data directory')
    # TODO: add a --log-level flag

    args = parser.parse_args(sys.argv[1:])

    run_id = os.path.basename(args.input_run_folder)
    fastq_dir = os.path.join(config['fastq_dir'], run_id)
    job_dir = os.path.join(config['jobs_dir'], run_id)
    
    logging.config.dictConfig(config.logging_config())
    main_logger = util.NamedAppLogger('main')

    setup_working_dirs(fastq_dir, job_dir)
    main_logger.info('Reading bcl data from ' + args.input_run_folder)
    main_logger.info('Fastq path is ' + fastq_dir)

    # Read RunInfo.xml and SampleSheet.csv for barcodes and masks
    main_logger.info('Reading the reads info from ' + args.input_run_folder)

    sample_sheet = analysis_driver.reader.sample_sheet.SampleSheet(args.input_run_folder)
    run_info = analysis_driver.reader.run_info.RunInfo(args.input_run_folder)
    if run_info.barcode_len:
        if not sample_sheet.check_barcodes() == run_info.barcode_len:
            main_logger.warn('Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' % (
                sample_sheet.check_barcodes(), run_info.barcode_len
                )
            )
    else:
        main_logger.warn('No barcode in RunInfo.xml')

    mask = run_info.mask.tostring(sample_sheet.check_barcodes())
    main_logger.info('bcl2fastq mask: ' + mask)  # example_mask = 'y150n,i6,y150n'

    if config['job_execution'] == 'pbs':
        run_pbs(
            logger=main_logger, input_run_folder=args.input_run_folder, job_dir=job_dir,
            run_id=run_id, fastq_dir=fastq_dir, mask=mask, sample_sheet=sample_sheet
        )

    util.demultiplex_feedback(run_id)

    main_logger.info('Done')
    return 0


def setup_working_dirs(*args):
    """
    Check for existence of working dirs and create them if necessary
    :param str args: Directories to check/create
    """
    for wd in args:
        if not os.path.exists(wd):
            logging_helper.debug('Creating: ' + wd)
            os.makedirs(wd)
        else:
            logging_helper.debug('Already exists: ' + wd)


def run_pbs(logger=None, input_run_folder=None, job_dir=None,
            run_id=None, fastq_dir=None, mask=None, sample_sheet=None):
    """
    Run PBS submission for compute jobs.
    :param logging.Logger logger: A main logger to use
    :param str input_run_folder: input_run_folder as passed to main
    :param str job_dir: The job dir for PBS scripts, BCBio, etc.
    :param str run_id: The basename of input_run_folder
    :param str fastq_dir: Expected path to generated fastq files
    :param str mask: A mask for bcl2fastq to use
    :param analysis_driver.reader.SampleSheet sample_sheet: The run's SampleSheet object
    """
    
    # Write bcl2fastq PBS script
    logger.info('Create bcl2fastq PBS script')
    bcl2fastq_pbs_name = os.path.join(job_dir, 'bcl2fastq_' + run_id + '.pbs')
    logger.info('bcl2fastq PBS file is ' + bcl2fastq_pbs_name)
    bcl2fastq_writer = writer.pbs_writer.BCL2FastqWriter(
        bcl2fastq_pbs_name, 'bcl2fastq', os.path.join(job_dir, 'bcl2fastq_pbs.log')
    )
    bcl2fastq_writer.write(job_dir, mask, input_run_folder, fastq_dir)

    # submit the bcl2fastq script to batch scheduler
    logger.info('Submitting: ' + os.path.join(job_dir, bcl2fastq_pbs_name))
    
    bcl2fastq_jobid = str(
        qsub_dependents.qsub([os.path.join(job_dir, bcl2fastq_pbs_name)])
    ).lstrip('b\'').rstrip('\'')
    logger.info('BCL2FASTQ jobId: ' + bcl2fastq_jobid)
    
    # Wait for fastqs to be created
    bcl2fastq_complete = os.path.join(job_dir, '.bcl2fastq_complete')
    logger.info('Waiting for creation of ' + bcl2fastq_complete)
    while not os.path.exists(bcl2fastq_complete):
        sleep(15)
    logger.info('bcl2fastq complete, executing fastqc and BCBio')
    
    # Write fastqc PBS script
    logger.info('Creating fastqc PBS script')
    fastqc_pbs_name = os.path.join(job_dir, 'fastqc_' + run_id + '.pbs')
    logger.info('Fastqc PBS file is ' + fastqc_pbs_name)

    sample_projects = list(sample_sheet.sample_projects.keys())
    fastqs = util.fastq_handler.flatten_fastqs(fastq_dir, sample_projects)
    fastqc_writer = writer.pbs_writer.FastqcWriter(
        fastqc_pbs_name, 'fastqc_' + run_id, os.path.join(job_dir, 'fastqc_pbs.log'), fastqs
    )
    fastqc_writer.write()
    
    # submit the fastqc script to the batch scheduler
    logger.info('Submitting: ' + os.path.join(job_dir, fastqc_pbs_name))
    fastqc_jobid = str(
        qsub_dependents.qsub([os.path.join(job_dir, fastqc_pbs_name)])
    ).lstrip('b\'').rstrip('\'')
    logger.info('FASTQC jobId: ' + fastqc_jobid)
    
    csv_writer = writer.BCBioSamplePrep(fastq_dir, job_dir, sample_sheet)
    samples = csv_writer.write()

    os.chdir(job_dir)
    bcbio_pbs = os.path.join(job_dir, 'bcbio_jobs.pbs')
    bcbio_writer = writer.pbs_writer.BCBioWriter(
        bcbio_pbs,
        'bcbio_' + run_id,
        'bcbio.log',
        len(samples)
    )

    for sample_id, csv_file, fastqs in samples:
        logger.info('Preparing samples for ' + sample_id)
        merged_csv = util.bcbio_prepare_samples(
            os.path.join(os.path.dirname(config['bcbio']), 'bcbio_prepare_samples.py'),
            csv_file
        )
        logger.info('Merged csv: ' + merged_csv)

        logger.info('Setting up BCBio run for ' + sample_id)
        util.setup_bcbio_run(
            config['bcbio'],
            os.path.join(config['location'], 'etc', 'bcbio_alignment.yaml'),
            os.path.join(job_dir, 'bcbio'),
            csv_file,
            *fastqs
        )

        bcbio_writer.add_bcbio_job(
            run_yaml=os.path.join(job_dir, 'samples_' + sample_id, 'config', 'samples_' + sample_id + '.yaml'),
            workdir=os.path.join(job_dir, 'samples_' + sample_id, 'work')
        )

    bcbio_writer.write()
    
    bcbio_jobid = str(
        # qsub_dependents.qsub_dependents([bcbio_pbs], jobid=fastqc_jobid)
        qsub_dependents.qsub([bcbio_pbs])
    ).lstrip('b\'').rstrip('\'')
    logger.info('BCBio jobId: ' + bcbio_jobid)
