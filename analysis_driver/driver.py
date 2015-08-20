import logging
import os
from analysis_driver import reader, writer, util, executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg  # imports the default config singleton

app_logger = logging.getLogger('driver')


def pipeline(input_run_folder):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
    """
    run_id = os.path.basename(input_run_folder)
    fastq_dir = os.path.join(cfg['fastq_dir'], run_id)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)

    app_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    app_logger.info('Fastq dir: ' + fastq_dir)
    app_logger.info('Job dir: ' + job_dir)

    setup_working_dirs(fastq_dir, job_dir)

    sample_sheet = reader.SampleSheet(input_run_folder)
    sample_sheet.validate()

    mask = sample_sheet.generate_mask()
    app_logger.info('bcl2fastq mask: ' + mask)  # example_mask = 'y150n,i6,y150n'

    if cfg['job_execution'] == 'pbs':
        run_pbs(
            logger=app_logger, input_run_folder=input_run_folder, job_dir=job_dir,
            run_id=run_id, fastq_dir=fastq_dir, mask=mask, sample_sheet=sample_sheet
        )

    app_logger.info('Done')
    return 0


def setup_working_dirs(*args):
    """
    Check for existence of working dirs and create them if necessary
    :param str args: Directories to check/create
    """
    for wd in args:
        if not os.path.exists(wd):
            app_logger.debug('Creating: ' + wd)
            os.makedirs(wd)
        else:
            app_logger.debug('Already exists: ' + wd)


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
    bcl2fastq_writer = writer.PBSWriter(
        bcl2fastq_pbs_name, 24, 12, 32, 'bcl2fastq_' + run_id, os.path.join(job_dir, 'bcl2fastq.log')
    )
    bcl2fastq_writer.write_line(writer.command_writer.bcl2fastq(mask, input_run_folder, fastq_dir))
    bcl2fastq_writer.save()

    # submit the bcl2fastq script to batch scheduler
    logger.info('Submitting ' + bcl2fastq_pbs_name)
    bcl2fastq_executor = executor.ClusterExecutor(os.path.join(job_dir, bcl2fastq_pbs_name), block=True)
    bcl2fastq_executor.start()

    bcl2fastq_exit_status = bcl2fastq_executor.join()
    logger.info('Exit status: ' + str(bcl2fastq_exit_status))
    if bcl2fastq_exit_status:
        raise AnalysisDriverError('Bcl2fastq failed')

    # Write fastqc PBS script
    fastqc_pbs_name = os.path.join(job_dir, 'fastqc_' + run_id + '.pbs')
    logger.info('Writing fastqc PBS script ' + fastqc_pbs_name)

    sample_projects = list(sample_sheet.sample_projects.keys())
    fastqs = util.fastq_handler.flatten_fastqs(fastq_dir, sample_projects)

    fastqc_writer = writer.PBSWriter(
        fastqc_pbs_name, 6, 8, 3, 'fastqc', os.path.join(job_dir, 'fastqc.log'), array=len(fastqs)
    )
    fastqc_writer.start_array()
    for idx, fastq in enumerate(fastqs):
        fastqc_writer.write_array_cmd(idx + 1, writer.command_writer.fastqc(fastq))
    fastqc_writer.finish_array()
    fastqc_writer.save()

    # TODO: writer should decide itself what the execution type is, and whether it needs to write a job array

    # submit the fastqc script to the batch scheduler
    logger.info('Submitting: ' + os.path.join(job_dir, fastqc_pbs_name))
    fastqc_executor = executor.ClusterExecutor(os.path.join(job_dir, fastqc_pbs_name))
    fastqc_executor.start()

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
            bcbio_csv_file = writer.write_bcbio_csv(job_dir, sample_id, id_fastqs)

            util.bcbio_prepare_samples(bcbio_csv_file)
            util.setup_bcbio_run(
                os.path.join(os.path.dirname(__file__), '..', 'etc', 'bcbio_alignment.yaml'),
                os.path.join(job_dir, 'bcbio'),
                bcbio_csv_file,
                *fastqs
            )

    bcbio_writer = writer.PBSWriter(
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

    bcbio_executor = executor.ClusterExecutor(bcbio_pbs)
    bcbio_executor.start()
