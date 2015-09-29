import os
from analysis_driver import reader, writer, util, executor
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg
from analysis_driver.notification import default as ntf

app_logger = get_logger('driver')


def pipeline(input_run_folder):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    """
    exit_status = 0
    run_id = os.path.basename(input_run_folder)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)
    fastq_dir = os.path.join(job_dir, 'fastq')
    app_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    app_logger.info('Fastq dir: ' + fastq_dir)
    app_logger.info('Job dir: ' + job_dir)

    run_info = reader.RunInfo(input_run_folder)
    samplesheet_csv = os.path.join(input_run_folder, 'SampleSheet.csv')
    if not run_info.mask.barcode_len or not os.path.exists(samplesheet_csv):
        app_logger.info('No sample sheet or barcodes found. Running in phiX mode')
        return pipeline_phix(input_run_folder)

    ntf.start_stage('setup')
    reader.transform_sample_sheet(input_run_folder)
    sample_sheet = reader.SampleSheet(input_run_folder)
    if not sample_sheet.validate(run_info.mask):
        raise AnalysisDriverError('Validation failed. Check SampleSheet.csv and RunInfo.xml.')

    mask = sample_sheet.generate_mask(run_info.mask)
    app_logger.info('bcl2fastq mask: ' + mask)  # example_mask = 'y150n,i6,y150n'
    ntf.end_stage('setup')
    
    # bcl2fastq
    ntf.start_stage('bcl2fastq')
    exit_status += _run_bcl2fastq(input_run_folder, run_id, fastq_dir, samplesheet_csv, mask).join()
    ntf.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status
    
    # start fastqc and bcbio
    ntf.start_stage('fastqc')
    fastqc_executor = _run_fastqc(run_id, fastq_dir)
    ntf.start_stage('bcbio')
    bcbio_executor = _run_bcbio(run_id, fastq_dir, job_dir, sample_sheet)

    # wait for fastqc and bcbio to finish
    fastqc_exit_status = fastqc_executor.join()
    ntf.end_stage('fastqc', fastqc_exit_status)
    bcbio_exit_status = bcbio_executor.join()
    ntf.end_stage('bcbio', bcbio_exit_status)

    # sort out exit statuses
    if bcbio_exit_status:
        return bcbio_exit_status
    exit_status += fastqc_exit_status + bcbio_exit_status
    
    # transfer output data
    ntf.start_stage('data_transfer')
    transfer_exit_status = _output_data(sample_sheet, job_dir)
    ntf.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status

    return exit_status


def pipeline_phix(input_run_folder):
    exit_status = 0
    run_id = os.path.basename(input_run_folder)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)
    fastq_dir = os.path.join(job_dir, 'fastq')

    # bcl2fastq
    ntf.start_stage('bcl2fastq')
    exit_status += _run_bcl2fastq(input_run_folder, run_id, fastq_dir).join()
    ntf.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status

    # fastqc
    ntf.start_stage('fastqc')
    exit_status += _run_fastqc(run_id, fastq_dir).join()
    ntf.end_stage('fastqc', exit_status)

    return exit_status


def _run_bcl2fastq(input_run_folder, run_id, fastq_dir, sample_sheet=None, mask=None):
    bcl2fastq_writer = writer.get_script_writer(
        'bcl2fastq',
        run_id,
        walltime=32,
        cpus=8,
        mem=32
    )
    bcl2fastq_script = writer.write_jobs(
        bcl2fastq_writer,
        [writer.bash_commands.bcl2fastq(input_run_folder, fastq_dir, sample_sheet, mask)]
    )

    app_logger.info('Submitting ' + bcl2fastq_script)
    bcl2fastq_executor = executor.ClusterExecutor(bcl2fastq_script, block=True)
    bcl2fastq_executor.start()
    return bcl2fastq_executor


def _run_fastqc(run_id, fastq_dir):
    fastqs = util.fastq_handler.find_all_fastqs(fastq_dir)
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
        [writer.bash_commands.fastqc(fq) for fq in fastqs]
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
        for sample_id, id_obj in proj_obj.sample_ids.items():
            id_fastqs = util.fastq_handler.find_fastqs(fastq_dir, sample_project, sample_id)
            merged_fastqs = util.bcbio_prepare_samples(job_dir, sample_id, id_fastqs)
            if merged_fastqs: 
                util.setup_bcbio_run(
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), '..', 'etc', 'bcbio_alignment.yaml'
                    ),
                    os.path.join(job_dir, 'bcbio'),
                    os.path.join(job_dir, 'samples_' + sample_id + '-merged.csv'),
                    merged_fastqs
                )
                bcbio_array_cmds.append(
                    writer.bash_commands.bcbio(
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
                        ),
                        threads=10
                    )
                )
            else:
                app_logger.warning('fastq merge failed - are some files missing or empty?')

    bcbio_writer = writer.get_script_writer(
        'bcbio',
        run_id,
        walltime=96,
        cpus=10,
        mem=64,
        jobs=len(bcbio_array_cmds)
    )
    for cmd in writer.bash_commands.bcbio_env_vars():
        bcbio_writer.write_line(cmd)

    bcbio_script = writer.write_jobs(
        bcbio_writer,
        bcbio_array_cmds,
        log_file_base=os.path.join(job_dir, 'bcbio')
    )

    bcbio_executor = executor.ClusterExecutor(bcbio_script, block=True)
    bcbio_executor.start()

    os.chdir(original_dir)

    return bcbio_executor


def _output_data(sample_sheet, job_dir):
    exit_status = 0
    for name, sample_project in sample_sheet.sample_projects.items():
        for sample_id in sample_project.sample_ids:

            output_dir = os.path.join(cfg['output_dir'], name, sample_id)
            try:
                os.makedirs(output_dir)
            except FileExistsError:
                pass

            bcbio_source_dir = os.path.join(job_dir, 'samples_' + sample_id + '-merged', 'final', sample_id)
            merged_fastq_dir = os.path.join(job_dir, 'merged')
            source_path_mapping = {
                'vcf': bcbio_source_dir,
                'bam': bcbio_source_dir,
                'fastq': merged_fastq_dir
            }

            ex = util.transfer_output_files(sample_id, output_dir, source_path_mapping)
            exit_status += ex

    return exit_status
