import os
from glob import glob
from analysis_driver import reader, writer, util, executor, clarity
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg  # imports the default config singleton
from analysis_driver.quality_control.genotype_validation import GenotypeValidation
from analysis_driver.notification import default as ntf

app_logger = get_logger('driver')


def pipeline(input_run_folder):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
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
    exit_status += executor.execute(
        [writer.bash_commands.bcl2fastq(input_run_folder, fastq_dir, sample_sheet.filename, mask)],
        job_name='bcl2fastq',
        run_id=run_id,
        walltime=32,
        cpus=8,
        mem=32
    ).join()
    # exit_status += _run_bcl2fastq(input_run_folder, run_id, fastq_dir, samplesheet_csv, mask).join()
    ntf.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status
    
    # fastqc
    ntf.start_stage('fastqc')
    # fastqc_executor = _run_fastqc(run_id, fastq_dir)
    fastqc_executor = executor.execute(
        [writer.bash_commands.fastqc(fq) for fq in util.fastq_handler.find_all_fastqs(fastq_dir)],
        job_name='fastqc',
        run_id=run_id,
        walltime=6,
        cpus=4,
        mem=2
    )

    valid_lanes = clarity.get_valid_lanes_from_HiseqX(run_info.flowcell_name)

    # merge fastq files
    ntf.start_stage('merge fastqs')
    sample_to_fastq_files = _bcio_prepare_sample(fastq_dir, job_dir, sample_sheet, valid_lanes)
    ntf.end_stage('merge fastqs')

    # genotype validation
    ntf.start_stage('genotype validation')
    genotype_validation = GenotypeValidation(sample_to_fastq_files, run_id)
    genotype_validation.start()

    # bcbio
    ntf.start_stage('bcbio')
    bcbio_executor = _run_bcbio(run_id, job_dir, sample_to_fastq_files)

    # wait for genotype_validation fastqc and bcbio to finish
    genotype_results = genotype_validation.join()
    app_logger.info('Written files: ' + str(genotype_results))
    ntf.end_stage('genotype validation')

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
    transfer_exit_status = _output_data(sample_sheet, job_dir, cfg['output_dir'], cfg['output_files'])
    ntf.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status

    return exit_status


def pipeline_phix(input_run_folder):
    # TODO: remove this, find a better way of having alternate pipelines
    exit_status = 0
    run_id = os.path.basename(input_run_folder)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)
    fastq_dir = os.path.join(job_dir, 'fastq')

    # bcl2fastq
    ntf.start_stage('bcl2fastq')
    exit_status += executor.execute(
        [writer.bash_commands.bcl2fastq(input_run_folder, fastq_dir)],
        job_name='bcl2fastq',
        run_id=run_id,
        walltime=32,
        cpus=8,
        mem=32
    ).join()
    ntf.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status

    # fastqc
    ntf.start_stage('fastqc')
    exit_status += executor.execute(
        [writer.bash_commands.fastqc(fq) for fq in util.fastq_handler.find_all_fastqs(fastq_dir)],
        job_name='fastqc',
        run_id=run_id,
        walltime=6,
        cpus=4,
        mem=2
    ).join()
    ntf.end_stage('fastqc', exit_status)

    return exit_status


def _bcio_prepare_sample(fastq_dir, job_dir, sample_sheet, valid_lanes=None):
    """
    Merge the fastq files per sample using bcbio prepare sample
    """
    sample_name_to_fastqs = {}
    for sample_project, proj_obj in sample_sheet.sample_projects.items():

        for sample_id, id_obj in proj_obj.sample_ids.items():
            if valid_lanes:
                #Only retrive the lanes that are considered valid
                fastq_files = []
                for lane in valid_lanes:
                    fastq_files.extend(util.fastq_handler.find_fastqs(fastq_dir, sample_project, sample_id, lane))
            else:
                #No information about valid lanes: don't filter anything
                fastq_files.extend(util.fastq_handler.find_fastqs(fastq_dir, sample_project, sample_id))
            merged_fastqs = util.bcbio_prepare_samples(
                job_dir,
                sample_id,
                util.fastq_handler.find_fastqs(fastq_dir, sample_project, sample_id)
            )
            sample_name_to_fastqs[sample_id] = merged_fastqs
    return sample_name_to_fastqs


def _run_bcbio(run_id, job_dir, sample_name_to_fastqs):
    run_template = os.path.join(
        os.path.dirname(__file__),
        '..', 'etc', 'bcbio_alignment_' + cfg['genome'] + '.yaml'
    )
    if not os.path.isfile(run_template):
        raise AnalysisDriverError(
            'Could not find BCBio run template ' + run_template + '. Is the correct genome set?'
        )

    original_dir = os.getcwd()
    os.chdir(job_dir)

    sample_preps = []
    bcbio_array_cmds = []
    for sample_id in sample_name_to_fastqs:
        sample_prep = [
            os.path.join(cfg['bcbio'], 'bin', 'bcbio_nextgen.py'),
            '-w template',
            run_template,
            os.path.join(job_dir, 'samples_' + sample_id + '-merged'),
            os.path.join(job_dir, 'samples_' + sample_id + '-merged.csv')
        ] + sample_name_to_fastqs.get(sample_id)
        sample_preps.append(' '.join(sample_prep))

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
                threads=16
            )
        )
    prep_status = executor.execute(sample_preps, env='local').join()
    app_logger.info('BCBio sample prep exit status: ' + str(prep_status))

    bcbio_executor = executor.execute(
        bcbio_array_cmds,
        prelim_cmds=writer.bash_commands.bcbio_env_vars(),
        job_name='bcbio',
        run_id=run_id,
        walltime=96,
        cpus=8,
        mem=64
    )
    os.chdir(original_dir)

    return bcbio_executor


def _output_data(sample_sheet, job_dir, output_dir, output_config):
    exit_status = 0
    for name, sample_project in sample_sheet.sample_projects.items():
        for sample_id in sample_project.sample_ids:

            output_loc = os.path.join(output_dir, name, sample_id)
            if not os.path.isdir(output_loc):
                os.makedirs(output_loc)

            for output_record in output_config:
                src_pattern = os.path.join(
                    job_dir,
                    os.path.join(*output_record['location']),
                    output_record['basename']
                ).replace('%s', sample_id)

                sources = glob(src_pattern)
                if sources:
                    source = sources[-1]

                    dest = os.path.join(
                        output_loc,
                        output_record.get('new_name', os.path.basename(source))
                    ).replace('%s', sample_id)
                    exit_status += util.transfer_output_file(
                        source,
                        dest
                    )

                else:
                    app_logger.warning('No files found for pattern ' + src_pattern)
                    exit_status += 1

    return exit_status
