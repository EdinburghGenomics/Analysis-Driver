import os
import logging
import logging.config
import argparse
from time import sleep

from analysis_driver import writer, util
from analysis_driver.legacy import qsub_dependents
import analysis_driver.reader
import analysis_driver.config as cfg

logging_helper = logging.getLogger('driver')


def main(argv):
    """
    :args:
        input_run_folder - absolute path to the input data directory, e.g. /abs/path/to/INPUT_DATA/run_id
    :return: exit status
    :rtype: int
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('input_run_folder', type=str, help='An absolute path to an input data directory')

    args = parser.parse_args(argv)
    config = cfg.Configuration()

    run_id = os.path.basename(args.input_run_folder)
    fastq_dir = os.path.join(config['fastq_dir'], run_id)
    job_dir = os.path.join(config['jobs_dir'], run_id)

    logging.config.dictConfig(config.logging())
    main_logger = logging.getLogger('main')

    main_logger.debug('lol')

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
            logger=main_logger, config=config, input_run_folder=args.input_run_folder, job_dir=job_dir,
            run_id=run_id, fastq_dir=fastq_dir, mask=mask, sample_sheet=sample_sheet
        )

    main_logger.info('Done')
    return 0


def setup_working_dirs(*args):
    for wd in args:
        if not os.path.exists(wd):
            logging_helper.debug('Creating: ' + wd)
            os.makedirs(wd)
        else:
            logging_helper.debug('Already exists: ' + wd)


def run_pbs(logger=None, config=None, input_run_folder=None, job_dir=None,
            run_id=None, fastq_dir=None, mask=None, sample_sheet=None):

    # Write bcl2fastq PBS script
    logger.info('Create bcl2fastq PBS script')
    bcl2fastq_pbs_name = os.path.join(job_dir, 'bcl2fastq_' + run_id + '.pbs')
    logger.info('bcl2fastq PBS file is ' + bcl2fastq_pbs_name)
    bcl2fastq_writer = writer.pbs_writer.BCL2FastqWriter(
        bcl2fastq_pbs_name, 'bcl2fastq', os.path.join(job_dir, 'bcl2fastq_pbs.log')
    )
    bcl2fastq_writer.write(mask, input_run_folder, fastq_dir)

    # Write fastqc PBS script
    logger.info('Creating fastqc PBS script')
    fastqc_pbs_name = os.path.join(job_dir, 'fastqc_' + run_id + '.pbs')
    logger.info('Fastqc PBS file is ' + fastqc_pbs_name)

    fastqc_writer = writer.pbs_writer.FastqcWriter(
        fastqc_pbs_name, 'fastqc', os.path.join(job_dir, 'fastqc_pbs.log')
    )
    fastqc_writer.write(fastq_dir, job_dir)

    # submit the bcl2fastq script to batch scheduler
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

    csv_writer = writer.BCBioCSVWriter(fastq_dir, job_dir, sample_sheet)
    csv_writer.write()
    fastqs = []

    for sample_project in sample_sheet.sample_projects:
        fastqs = fastqs + util.find_fastqs(os.path.join(fastq_dir, sample_project))

    os.chdir(job_dir)

    bcbio_pbs = os.path.join(job_dir, 'bcbio.pbs')
    bcbio_writer = writer.pbs_writer.BCBioWriter(bcbio_pbs, 'bcbio', 'bcbio.log')
    util.setup_bcbio_run(
        config['bcbio'],
        os.path.join(os.path.dirname(__file__), 'etc', 'bcbio_alignment.yaml'),
        os.path.join(job_dir, 'bcbio'),
        os.path.join(job_dir, 'samples.csv'),
        fastqs
    )

    bcbio_writer.write(
        config['bcbio'],
        os.path.join(job_dir, 'samples', 'config', 'samples.yaml'),
        os.path.join(job_dir, 'samples', 'work')
    )
    bcbio_jobid = str(
        # qsub_dependents.qsub_dependents([bcbio_pbs], jobid=fastqc_jobid)
        qsub_dependents.qsub([bcbio_pbs])
    ).lstrip('b\'').rstrip('\'')
    logger.info('BCBio jobId: ' + bcbio_jobid)


if __name__ == '__main__':
    from sys import argv
    main(argv)
