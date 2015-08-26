import os
from analysis_driver import reader, writer, util, executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg  # imports the default config singleton

app_logger = get_logger('driver')


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

    reader.transform_sample_sheet(input_run_folder)

    sample_sheet = reader.SampleSheet(input_run_folder)
    sample_sheet.validate()

    mask = sample_sheet.generate_mask()
    app_logger.info('bcl2fastq mask: ' + mask)  # example_mask = 'y150n,i6,y150n'

    bcl2fastq_writer = writer.get_script_writer(
        'bcl2fastq',
        run_id,
        walltime=24,
        cpus=12,
        mem=32
    )
    bcl2fastq_script = writer.write_jobs(
        bcl2fastq_writer,
        [writer.commands.bcl2fastq(mask, input_run_folder, fastq_dir)]
    )

    app_logger.info('Submitting ' + bcl2fastq_script)
    bcl2fastq_executor = executor.ClusterExecutor(bcl2fastq_script, block=True)
    bcl2fastq_executor.start()
    bcl2fastq_exit_status = bcl2fastq_executor.join()
    app_logger.info('Exit status: ' + str(bcl2fastq_exit_status))
    if bcl2fastq_exit_status:
        raise AnalysisDriverError('Bcl2fastq failed')

    sample_projects = list(sample_sheet.sample_projects.keys())
    fastqs = util.fastq_handler.flatten_fastqs(fastq_dir, sample_projects)

    fastqc_writer = writer.get_script_writer(
        'fastqc',
        run_id,
        walltime=6,
        cpus=8,
        mem=3,
        jobs=len(fastqs)
    )
    fastqc_script = writer.write_jobs(
        fastqc_writer,
        [writer.commands.fastqc(fq) for fq in fastqs]
    )
    app_logger.info('Submitting: ' + fastqc_script)
    fastqc_executor = executor.ClusterExecutor(fastqc_script, block=True)
    fastqc_executor.start()

    os.chdir(job_dir)

    bcbio_array_cmds = []
    for sample_project, proj_obj in sample_sheet.sample_projects.items():
        proj_fastqs = util.fastq_handler.find_fastqs(fastq_dir, sample_project)

        for sample_id, id_obj in proj_obj.sample_ids.items():

            bcbio_array_cmds.append(
                writer.commands.bcbio(
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
                os.path.join(job_dir, 'samples_' + sample_id + '-merged.csv'),
                *id_fastqs
            )

    bcbio_writer = writer.get_script_writer(
        'bcbio',
        run_id,
        walltime=72,
        cpus=8,
        mem=64,
        jobs=len(bcbio_array_cmds)
    )
    for cmd in writer.commands.bcbio_java_paths():
        bcbio_writer.write_line(cmd)

    bcbio_script = writer.write_jobs(
        bcbio_writer,
        bcbio_array_cmds
    )

    bcbio_executor = executor.ClusterExecutor(bcbio_script, block=True)
    bcbio_executor.start()

    fastqc_exit_status = fastqc_executor.join()
    bcbio_exit_status = bcbio_executor.join()

    app_logger.info('rsync goes here')
    # util.transfer_output_data(os.path.basename(input_run_folder))

    app_logger.info('bcl2fastq exit status: ' + str(bcl2fastq_exit_status))
    app_logger.info('fastqc exit status: ' + str(fastqc_exit_status))
    app_logger.info('bcbio exit status: ' + str(bcbio_exit_status))
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
