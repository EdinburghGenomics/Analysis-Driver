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
    main_logger = logging.getLogger('main')

    run_id = os.path.basename(input_run_folder)
    fastq_dir = os.path.join(cfg['fastq_dir'], run_id)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)

    main_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    main_logger.info('Fastq dir: ' + fastq_dir)
    main_logger.info('Job dir: ' + job_dir)

    setup_working_dirs(fastq_dir, job_dir)

    sample_sheet = reader.sample_sheet.SampleSheet(input_run_folder)
    run_info = reader.run_info.RunInfo(input_run_folder)
    validate_run_info(run_info, sample_sheet)

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


def validate_run_info(run_info, sample_sheet):
    assert run_info.barcode_len, 'No barcode found in RunInfo.xml'
    assert sample_sheet.check_barcodes() == run_info.barcode_len, 'Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' %\
        (sample_sheet.check_barcodes(), run_info.barcode_len)


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
    bcl2fastq_pbs_name = os.path.join(job_dir, 'bcl2fastq_' + run_id + '.pbs')
    logger.info('Writing bcl2fastq PBS script ' + bcl2fastq_pbs_name)
    bcl2fastq_writer = writer.pbs_script_writer.PBSWriter(
        bcl2fastq_pbs_name, 24, 12, 32, 'bcl2fastq_' + run_id, os.path.join(job_dir, 'bcl2fastq.log')
    )
    bcl2fastq_writer.write_line(writer.command_writer.bcl2fastq(mask, input_run_folder, fastq_dir))
    bcl2fastq_writer.save()

    # submit the bcl2fastq script to batch scheduler
    logger.info('Submitting ' + bcl2fastq_pbs_name)

    bcl2fastq_jobid = str(
        qsub_dependents.qsub([os.path.join(job_dir, bcl2fastq_pbs_name)])
    ).lstrip('b\'').rstrip('\'')
    logger.info('bcl2fastq job id: ' + bcl2fastq_jobid)

    # Wait for fastqs to be created
    bcl2fastq_complete = os.path.join(job_dir, '.bcl2fastq_complete')
    logger.info('Waiting for creation of ' + bcl2fastq_complete)
    while not os.path.exists(bcl2fastq_complete):
        sleep(15)
    logger.info('Bcl2fastq complete, proceeding')

    # Write fastqc PBS script
    fastqc_pbs_name = os.path.join(job_dir, 'fastqc_' + run_id + '.pbs')
    logger.info('Writing fastqc PBS script ' + fastqc_pbs_name)

    sample_projects = list(sample_sheet.sample_projects.keys())
    print(sample_sheet.sample_projects['10015AT'].sample_ids)
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

    os.chdir(job_dir)
    bcbio_pbs = os.path.join(job_dir, 'bcbio_jobs.pbs')

    bcbio_array_cmds = []
    for sample_project, proj_obj in sample_sheet.sample_projects.items():
        proj_fastqs = util.fastq_handler.find_fastqs(fastq_dir, sample_project)

        for sample_id, id_obj in proj_obj.sample_ids.items():

            bcbio_array_cmds.append(
                writer.command_writer.bcbio(
                    os.path.join(job_dir, 'samples_' + sample_id, 'config', 'samples_' + sample_id + '.yaml'),
                    os.path.join(job_dir, 'samples_' + sample_id, 'work')
                )
            )

            id_fastqs = proj_fastqs[sample_id]

            csv_writer = writer.BCBioCSVWriter(job_dir, sample_id, id_fastqs)

            util.bcbio_prepare_samples(csv_writer.csv_file)
            util.setup_bcbio_run(
                os.path.join(os.path.dirname(__file__), '..', 'etc', 'bcbio_alignment.yaml'),
                os.path.join(job_dir, 'bcbio'),
                csv_writer.csv_file,
                *fastqs
            )

    bcbio_writer = writer.pbs_script_writer.PBSWriter(
        bcbio_pbs,
        72,
        8,
        64,
        'bcbio',
        'bcbio.log',
        array=len(bcbio_array_cmds)
    )

    for cmd in writer.command_writer.bcbio_java_paths():
        bcbio_writer.write_line(cmd)
    bcbio_writer.start_array()

    for idx, cmd in enumerate(bcbio_array_cmds):
        bcbio_writer.write_array_cmd(idx + 1, cmd)
        logger.info('Written command number ' + str(idx))

    bcbio_writer.finish_array()
    bcbio_writer.save()

    bcbio_jobid = str(
        # qsub_dependents.qsub_dependents([bcbio_pbs], jobid=fastqc_jobid)
        qsub_dependents.qsub([bcbio_pbs])
    ).lstrip('b\'').rstrip('\'')
    logger.info('BCBio jobId: ' + bcbio_jobid)
