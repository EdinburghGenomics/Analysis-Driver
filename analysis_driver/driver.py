import os
import yaml
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
    app_logger.info('bcl2fastq mask: ' + mask)  # e.g: mask = 'y150n,i6,y150n'
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
    ntf.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status
    
    # fastqc
    ntf.start_stage('fastqc')
    fastqc_executor = executor.execute(
        [writer.bash_commands.fastqc(fq) for fq in util.fastq_handler.find_all_fastqs(fastq_dir)],
        job_name='fastqc',
        run_id=run_id,
        walltime=6,
        cpus=4,
        mem=2
    )

    valid_lanes = clarity.get_valid_lanes_from_HiseqX(run_info.flowcell_name)
    app_logger.info('Valid lanes: ' + str(valid_lanes))

    # merge fastq files
    ntf.start_stage('merge fastqs')
    sample_to_fastq_files = _bcbio_prepare_samples(fastq_dir, job_dir, sample_sheet, valid_lanes)
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


def _bcbio_prepare_samples(fastq_dir, job_dir, sample_sheet, valid_lanes=None):
    """
    Merge the fastq files per sample using bcbio prepare sample
    """
    sample_name_to_fastqs = {}
    for sample_project, proj_obj in sample_sheet.sample_projects.items():

        for sample_id, id_obj in proj_obj.sample_ids.items():
            fastq_files = []
            if valid_lanes:
                # only retrive the lanes that are considered valid
                for lane in valid_lanes:
                    fastq_files.extend(util.fastq_handler.find_fastqs(fastq_dir, sample_project, sample_id, lane))
            else:
                # no information about valid lanes: don't filter anything
                fastq_files.extend(util.fastq_handler.find_fastqs(fastq_dir, sample_project, sample_id))
            #Query the LIMS to get the user sample id
            user_sample_id = clarity.get_user_sample_name(sample_id)
            merged_fastqs = util.bcbio_prepare_samples(
                job_dir,
                sample_id,
                fastq_files,
                user_sample_id=user_sample_id
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
    app_logger.debug(str(sample_name_to_fastqs))

    sample_preps = []
    bcbio_array_cmds = []
    for sample_name in sample_name_to_fastqs:
        app_logger.debug('Setting up sample: ' + sample_name)

        bcbio_dir = os.path.join(job_dir, 'samples_' + sample_name + '-merged')

        sample_prep = [
            os.path.join(cfg['bcbio'], 'bin', 'bcbio_nextgen.py'),
            '-w template',
            run_template,
            bcbio_dir,
            bcbio_dir + '.csv'
        ] + sample_name_to_fastqs.get(sample_name)

        run_yaml = os.path.join(bcbio_dir, 'config', 'samples_' + sample_name + '-merged.yaml')
        bcbio_array_cmds.append(
            writer.bash_commands.bcbio(
                run_yaml,
                os.path.join(bcbio_dir, 'work'),
                threads=16
            )
        )
        prep_status = executor.execute([' '.join(sample_prep)], env='local').join()
        app_logger.info('BCBio sample prep exit status: ' + str(prep_status))

        user_sample_id = clarity.get_user_sample_name(sample_name)
        if user_sample_id:
            app_logger.debug('Found user sample: ' + user_sample_id)

            with open(run_yaml, 'r') as i:
                run_config = yaml.load(i)
            run_config['fc_name'] = user_sample_id
            with open(run_yaml, 'w') as o:
                o.write(yaml.safe_dump(run_config, default_flow_style=False))

    bcbio_executor = executor.execute(
        bcbio_array_cmds,
        prelim_cmds=writer.bash_commands.bcbio_env_vars(),
        job_name='bcbio',
        run_id=run_id,
        walltime=120,
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

            user_sample_id = clarity.get_user_sample_name(sample_id)
            if not user_sample_id:
                user_sample_id = sample_id

            for output_record in output_config:
                src_pattern = os.path.join(
                    job_dir,
                    os.path.join(*output_record['location']),
                    output_record['basename']
                ).format(runfolder=sample_id, sample_id=user_sample_id)

                sources = glob(src_pattern)
                if sources:
                    source = sources[-1]

                    dest = os.path.join(
                        output_loc,
                        output_record.get('new_name', os.path.basename(source))
                    ).format(sample_id=user_sample_id)
                    exit_status += util.transfer_output_file(
                        source,
                        dest
                    )

                else:
                    app_logger.warning('No files found for pattern ' + src_pattern)
                    exit_status += 1

            with open(os.path.join(output_loc, 'run_config.yaml'), 'w') as f:
                f.write(cfg.report())

    return exit_status
