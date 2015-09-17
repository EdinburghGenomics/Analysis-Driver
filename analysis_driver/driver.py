import os
from analysis_driver import reader, writer, util, executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.notification import default_notification_center as ntf
from analysis_driver.config import default as cfg  # imports the default config singleton
from analysis_driver.quality_control import genotype_validation
from analysis_driver.quality_control.genotype_validation import GenotypeValidation

app_logger = get_logger('driver')


def pipeline(input_run_folder):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
    """
    run_id = os.path.basename(input_run_folder)

    ntf.start_pipeline(run_id)
    ntf.start_stage('setup')

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

    ntf.end_stage('setup', run_id)

    # bcl2fastq
    ntf.start_stage('bcl2fastq')
    bcl2fastq_exit_status = _run_bcl2fastq(input_run_folder, run_id, fastq_dir, mask).join()
    
    ntf.end_stage('bcl2fastq', run_id, bcl2fastq_exit_status, stop_on_error=True)

    # fastqc
    ntf.start_stage('fastqc')
    fastqc_executor = _run_fastqc(run_id, fastq_dir, sample_sheet)

    #Merge the fastq files
    ntf.start_stage('merge fastq')
    sample_to_fastq_files = _bcio_prepare_sample(fastq_dir, job_dir, sample_sheet)
    ntf.end_stage('merge fastq', run_id)

    GenotypeValidation()
    # bcbio
    ntf.start_stage('bcbio')
    bcbio_executor = _run_bcbio(run_id, job_dir, sample_to_fastq_files)

    # wait for fastqc and bcbio to finish
    fastqc_exit_status = fastqc_executor.join()
    ntf.end_stage('fastqc', run_id, fastqc_exit_status)

    bcbio_exit_status = bcbio_executor.join()
    ntf.end_stage('bcbio', run_id, bcbio_exit_status)

    # rsync goes here
    # util.transfer_output_data(os.path.basename(input_run_folder))
    ntf.close()
    ntf.end_pipeline(run_id)
    return 0


def _run_bcl2fastq(input_run_folder, run_id, fastq_dir, mask):
    bcl2fastq_writer = writer.get_script_writer(
        'bcl2fastq',
        run_id,
        walltime=32,
        cpus=8,
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
        cpus=4,
        mem=2,
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

def _bcio_prepare_sample(fastq_dir, job_dir, sample_sheet):
    """
    Merge the fastq files per sample using bcbio prepare sample
    """
    sample_name_to_fastqs = {}
    for sample_project, proj_obj in sample_sheet.sample_projects.items():
        proj_fastqs = util.fastq_handler.find_fastqs(fastq_dir, sample_project)

        for sample_id, id_obj in proj_obj.sample_ids.items():
            merged_fastqs = util.bcbio_prepare_samples(job_dir, sample_id, proj_fastqs[sample_id])
            sample_name_to_fastqs[sample_id]=merged_fastqs
    return sample_name_to_fastqs


def _run_bcbio(run_id, job_dir, sample_name_to_fastqs):
    original_dir = os.getcwd()
    os.chdir(job_dir)

    bcbio_array_cmds = []
    list_fastq_files = []
    sample_names = []
    for sample_id in sample_name_to_fastqs:

        util.setup_bcbio_run(
            os.path.join(os.path.dirname(__file__), '..', 'etc', 'bcbio_alignment.yaml'),
            os.path.join(job_dir, 'bcbio'),
            os.path.join(job_dir, 'samples_' + sample_id + '-merged.csv'),
            *sample_name_to_fastqs.get(sample_id)
        )
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

    bcbio_writer = writer.get_script_writer(
        'bcbio',
        run_id,
        walltime=72,
        cpus=16,
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
