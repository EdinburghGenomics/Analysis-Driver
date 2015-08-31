import os
from analysis_driver import reader, writer, util, executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.notification import default_notification_center as ntf
from analysis_driver.config import default as cfg  # imports the default config singleton

app_logger = get_logger('driver')


def pipeline(input_run_folder):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
    """
    ntf.start_step('setup')

    run_id = os.path.basename(input_run_folder)
    fastq_dir = os.path.join(cfg['fastq_dir'], run_id)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)

    app_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    app_logger.info('Fastq dir: ' + fastq_dir)
    app_logger.info('Job dir: ' + job_dir)

    reader.transform_sample_sheet(input_run_folder)

    sample_sheet = reader.SampleSheet(input_run_folder)
    sample_sheet_errors = sample_sheet.validate()
    if sample_sheet_errors:
        raise AnalysisDriverError('Sample sheet validation failed. See log messages.')

    mask = sample_sheet.generate_mask()
    app_logger.info('bcl2fastq mask: ' + mask)  # example_mask = 'y150n,i6,y150n'

    ntf.finish_step('setup')

    # bcl2fastq
    ntf.start_step('bcl2fastq')
    bcl2fastq_exit_status = _run_bcl2fastq(input_run_folder, run_id, fastq_dir, mask).join()

    if bcl2fastq_exit_status:
        ntf.fail_step('bcl2fastq')
        raise AnalysisDriverError('Bcl2fastq failed')

    ntf.finish_step('bcl2fastq')

    # fastqc
    ntf.start_step('fastqc')
    fastqc_executor = _run_fastqc(run_id, fastq_dir, sample_sheet)

    # bcbio
    ntf.start_step('bcbio')
    bcbio_executor = _run_bcbio(run_id, fastq_dir, job_dir, sample_sheet)

    # wait for fastqc and bcbio to finish
    fastqc_exit_status = fastqc_executor.join()
    if fastqc_exit_status:
        ntf.fail_step('fastqc')
    else:
        ntf.finish_step('fastqc')

    bcbio_exit_status = bcbio_executor.join()
    if bcbio_exit_status:
        ntf.fail_step('bcbio')
    else:
        ntf.finish_step('bcbio')

    app_logger.info('rsync goes here')
    # util.transfer_output_data(os.path.basename(input_run_folder))

    app_logger.info('Done')
    return 0


def _run_bcl2fastq(input_run_folder, run_id, fastq_dir, mask):
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

    return bcl2fastq_executor


def _run_fastqc(run_id, fastq_dir, sample_sheet):

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

    return fastqc_executor


def _run_bcbio(run_id, fastq_dir, job_dir, sample_sheet):
    original_dir = os.getcwd()
    os.chdir(job_dir)

    bcbio_array_cmds = []
    for sample_project, proj_obj in sample_sheet.sample_projects.items():
        proj_fastqs = util.fastq_handler.find_fastqs(fastq_dir, sample_project)

        for sample_id, id_obj in proj_obj.sample_ids.items():

            bcbio_array_cmds.append(
                writer.commands.bcbio(
                    os.path.join(
                        job_dir,
                        'samples_' + sample_id + '-merged',
                        'config',
                        'samples_' + sample_id + '-merged.yaml'
                    ),
                    os.path.join(
                        job_dir,
                        'samples_' + sample_id + '-merged',
                        'work'
                    )
                )
            )

            id_fastqs = proj_fastqs[sample_id]

            merged_fastqs = util.bcbio_prepare_samples(job_dir, sample_id, id_fastqs)
            #TODO: check that the merged fastqs exists
            util.setup_bcbio_run(
                os.path.join(os.path.dirname(__file__), '..', 'etc', 'bcbio_alignment.yaml'),
                os.path.join(job_dir, 'bcbio'),
                os.path.join(job_dir, 'samples_' + sample_id + '-merged.csv'),
                *merged_fastqs
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

    os.chdir(original_dir)

    return bcbio_executor
