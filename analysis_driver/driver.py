import logging
from time import sleep

import os

from analysis_driver import reader, writer, util
from analysis_driver.legacy import qsub_dependents
from analysis_driver.config import default as cfg  # imports the default config singleton

logging_helper = logging.getLogger('driver')


def pipeline(input_run_folder):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
    """

    run_id = os.path.basename(input_run_folder)
    fastq_dir = os.path.join(cfg['fastq_dir'], run_id)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)

    main_logger = logging.getLogger('main')

    setup_working_dirs(fastq_dir, job_dir)
    main_logger.info('Reading bcl data from ' + input_run_folder)
    main_logger.info('Fastq path is ' + fastq_dir)

    # Read RunInfo.xml and SampleSheet.csv for barcodes and masks
    main_logger.info('Reading the reads info from ' + input_run_folder)

    sample_sheet = reader.sample_sheet.SampleSheet(input_run_folder)
    run_info = reader.run_info.RunInfo(input_run_folder)
    if run_info.barcode_len:
        if not sample_sheet.check_barcodes() == run_info.barcode_len:
            main_logger.warn(
                'Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' % (
                    sample_sheet.check_barcodes(), run_info.barcode_len
                )
            )
    else:
        main_logger.warn('No barcode in RunInfo.xml')

    mask = run_info.mask.tostring(sample_sheet.check_barcodes())
    main_logger.info('bcl2fastq mask: ' + mask)  # example_mask = 'y150n,i6,y150n'

    if cfg['job_execution'] == 'pbs':
        run_pbs(
            logger=main_logger, input_run_folder=input_run_folder, job_dir=job_dir,
            run_id=run_id, fastq_dir=fastq_dir, mask=mask, sample_sheet=sample_sheet
        )

    # util.demultiplex_feedback(run_id)

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
    bcl2fastq_writer = writer.pbs_script_writer.PBSWriter(
        bcl2fastq_pbs_name, 24, 12, 32, 'bcl2fastq', os.path.join(job_dir, 'bcl2fastq.log')
    )
    bcl2fastq_writer.write_line(writer.command_writer.bcl2fastq(mask, input_run_folder, fastq_dir))
    bcl2fastq_writer.save()

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

    fastqc_writer = writer.pbs_script_writer.PBSWriter(
        fastqc_pbs_name, 6, 8, 3, 'fastqc', os.path.join(job_dir, 'fastqc.log'), array=len(fastqs)
    )
    for idx, fastq in enumerate(fastqs):
        fastqc_writer.write_array_cmd(idx + 1, writer.command_writer.fastqc(fastq))
    fastqc_writer.finish_array()
    fastqc_writer.save()
    
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

    bcbio_writer = writer.pbs_script_writer.PBSWriter(
        bcbio_pbs, 72, 8, 64, 'bcbio', 'bcbio.log', array=len(samples)
    )

    for idx, (sample_id, csv_file, fastqs) in enumerate(samples):
        logger.info('Preparing samples for ' + sample_id)
        merged_csv = util.bcbio_prepare_samples(
            os.path.join(os.path.dirname(cfg['bcbio']), 'bcbio_prepare_samples.py'),
            csv_file
        )
        logger.info('Merged csv: ' + merged_csv)

        logger.info('Setting up BCBio run for ' + sample_id)
        util.setup_bcbio_run(
            cfg['bcbio'],
            os.path.join(os.path.dirname(os.path.dirname(__file__)), 'etc', 'bcbio_alignment.yaml'),
            os.path.join(job_dir, 'bcbio'),
            csv_file,
            *fastqs
        )

        bcbio_writer.write_array_cmd(
            idx + 1,
            writer.command_writer.bcbio(
                os.path.join(job_dir, 'samples_' + sample_id, 'config', 'samples_' + sample_id + '.yaml'),
                os.path.join(job_dir, 'samples_' + sample_id, 'work')
            )
        )

    bcbio_writer.finish_array()
    bcbio_writer.save()
    
    bcbio_jobid = str(
        # qsub_dependents.qsub_dependents([bcbio_pbs], jobid=fastqc_jobid)
        qsub_dependents.qsub([bcbio_pbs])
    ).lstrip('b\'').rstrip('\'')
    logger.info('BCBio jobId: ' + bcbio_jobid)
