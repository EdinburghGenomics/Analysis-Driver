import os
import time
import yaml
from glob import glob
import shutil
from analysis_driver import reader, writer, util, executor, clarity
from analysis_driver.dataset_scanner import RunDataset, SampleDataset
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg  # imports the default config singleton
from analysis_driver.notification import default as ntf
from analysis_driver.report_generation.run_report import RunCrawler
from analysis_driver.transfer_data import prepare_run_data, prepare_sample_data, output_sample_data, output_run_data

app_logger = get_logger('driver')

def pipeline(dataset):

    if isinstance(dataset, RunDataset):
        exit_status = demultiplexing_pipeline(dataset)
    elif isinstance(dataset, SampleDataset):
        exit_status = variant_calling_pipeline(dataset)

    if exit_status != 0:
        dataset.fail()
    else:
        dataset.succeed()
    return exit_status

def demultiplexing_pipeline(dataset):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
    """
    exit_status = 0
    ntf.start_stage('transfer')
    input_run_folder = prepare_run_data(dataset)
    ntf.end_stage('transfer')

    dataset.start()

    run_id = os.path.basename(input_run_folder)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)
    fastq_dir = os.path.join(job_dir, 'fastq')
    app_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    app_logger.info('Fastq dir: ' + fastq_dir)
    app_logger.info('Job dir: ' + job_dir)

    run_info = reader.RunInfo(input_run_folder)

    ntf.start_stage('setup')
    reader.transform_sample_sheet(input_run_folder)
    sample_sheet = reader.SampleSheet(input_run_folder)
    if not sample_sheet.validate(run_info.mask):
        raise AnalysisDriverError('Validation failed. Check SampleSheet.csv and RunInfo.xml.')

    mask = sample_sheet.generate_mask(run_info.mask)
    app_logger.info('bcl2fastq mask: ' + mask)  # e.g: mask = 'y150n,i6,y150n'
    #Send the information about the run to the rest API
    crawler = RunCrawler(run_id, sample_sheet)
    crawler.send_data()
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
        cpus=1,
        mem=2
    )

    # md5sum
    ntf.start_stage('md5sum')
    md5sum_executor = executor.execute(
        [writer.bash_commands.md5sum(fq) for fq in util.fastq_handler.find_all_fastqs(fastq_dir)],
        job_name='md5sum',
        run_id=run_id,
        walltime=6,
        cpus=1,
        mem=2,
        log_command = False
    )

    valid_lanes = clarity.get_valid_lanes(run_info.flowcell_name)

    fastqc_exit_status = fastqc_executor.join()
    ntf.end_stage('fastqc', fastqc_exit_status)
    md5_exit_status = md5sum_executor.join()
    ntf.end_stage('md5sum', md5_exit_status)

    #copy the Samplesheet Runinfo.xml run_parameters.xml to the fastq dir
    for f in ['SampleSheet.csv', 'SampleSheet_analysis_driver.csv', 'runParameters.xml',
              'RunInfo.xml', 'RTAConfiguration.xml']:
        shutil.copy2(os.path.join(input_run_folder,f), os.path.join(fastq_dir, f))
    shutil.copytree(os.path.join(input_run_folder,'InterOp'), os.path.join(fastq_dir,'InterOp'))

    #Find conversion xml file and send the results to the rest API
    conversion_xml = os.path.join(fastq_dir, 'Stats','ConversionStats.xml')
    if os.path.exists(conversion_xml):
        crawler = RunCrawler(run_id, sample_sheet, conversion_xml)
        json_file = os.path.join(fastq_dir, 'demultiplexing_results.json')
        crawler.write_json(json_file)
        sample_dir = os.path.join(cfg['output_dir'],'samples')
        os.makedirs(sample_dir, exist_ok=True)
        crawler.update_json_per_sample(sample_dir)
        crawler.send_data()
    else:
        app_logger.error('File not found: %s'%conversion_xml)
        exit_status+=1

    ntf.start_stage('data_transfer')
    transfer_exit_status = output_run_data(fastq_dir, run_id)
    ntf.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status + fastqc_exit_status + md5_exit_status

    if exit_status == 0:
        ntf.start_stage('cleanup')
        exit_status += _cleanup(run_id)
        ntf.end_stage('cleanup', exit_status)
    return exit_status



def variant_calling_pipeline(dataset):
    """
    :param str input_run_folder: Full path to an input data directory
    :return: Exit status
    :rtype: int
    """
    exit_status = 0
    dataset.start()
    fastq_files = prepare_sample_data(dataset)


    sample_id = dataset.name
    sample_dir = os.path.join(cfg['jobs_dir'], sample_id)
    app_logger.info('Job dir: ' + sample_dir)

    # merge fastq files
    ntf.start_stage('merge fastqs')
    fastq_pair = _bcbio_prepare_sample(sample_dir, sample_id, fastq_files)
    app_logger.debug('sample fastq files: ' + str(fastq_pair))
    ntf.end_stage('merge fastqs')

    # fastqc2
    ntf.start_stage('sample_fastqc')
    fastqc2_executor = executor.execute(
        [writer.bash_commands.fastqc(fastq_file) for fastq_file in fastq_pair],
        job_name='fastqc2',
        run_id=sample_id,
        walltime=10,
        cpus=1,
        mem=2
    )

    # genotype validation
    # ntf.start_stage('genotype validation')
    # genotype_validation = GenotypeValidation(sample_to_fastq_files, run_id)
    # genotype_validation.start()

    # bcbio
    ntf.start_stage('bcbio')
    bcbio_executor = _run_bcbio(sample_id, sample_dir, fastq_pair)

    # wait for genotype_validation fastqc and bcbio to finish
    # genotype_results = genotype_validation.join()
    # app_logger.info('Written files: ' + str(genotype_results))
    # ntf.end_stage('genotype validation')

    fastqc2_exit_status = fastqc2_executor.join()
    ntf.end_stage('sample_fastqc2', fastqc2_exit_status)

    bcbio_exit_status = bcbio_executor.join()
    ntf.end_stage('bcbio', bcbio_exit_status)

    # sort out exit statuses
    if bcbio_exit_status:
        return bcbio_exit_status
    exit_status += fastqc2_exit_status + bcbio_exit_status

    # transfer output data
    ntf.start_stage('data_transfer')
    transfer_exit_status = output_sample_data(sample_id, sample_dir, cfg['output_dir'], cfg['output_files'])
    ntf.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status

    if exit_status == 0:
        ntf.start_stage('cleanup')
        exit_status += _cleanup(sample_id)
        ntf.end_stage('cleanup', exit_status)
    return exit_status


def _bcbio_prepare_sample(job_dir, sample_id, fastq_files):
    """
    Merge the fastq files per sample using bcbio prepare sample
    """
    user_sample_id = clarity.get_user_sample_name(sample_id)
    if not user_sample_id:
        user_sample_id = sample_id
    cmd = util.bcbio_prepare_samples_cmd(
            job_dir,
            sample_id,
            fastq_files,
            user_sample_id=user_sample_id
          )

    exit_status = executor.execute([cmd], job_name='bcbio_prepare_samples', run_id=sample_id).join()
    sample_fastqs = glob(os.path.join(job_dir, 'merged', user_sample_id + '_R?.fastq.gz'))

    app_logger.info('bcbio_prepare_samples finished with exit status ' + str(exit_status))

    return sample_fastqs

def _run_bcbio(sample_id, sample_dir, sample_fastqs):
    run_template = os.path.join(
        os.path.dirname(__file__),
        '..', 'etc', 'bcbio_alignment_' + cfg['genome'] + '.yaml'
    )
    if not os.path.isfile(run_template):
        raise AnalysisDriverError(
            'Could not find BCBio run template ' + run_template + '. Is the correct genome set?'
        )

    original_dir = os.getcwd()
    os.chdir(sample_dir)
    app_logger.debug(str(sample_fastqs))

    app_logger.debug('Setting up sample: ' + sample_id)

    bcbio_dir = os.path.join(sample_dir, 'samples_' + sample_id + '-merged')

    sample_prep = [
        os.path.join(cfg['bcbio'], 'bin', 'bcbio_nextgen.py'),
        '-w template',
        run_template,
        bcbio_dir,
        bcbio_dir + '.csv'
    ] + sample_fastqs

    run_yaml = os.path.join(bcbio_dir, 'config', 'samples_' + sample_id + '-merged.yaml')
    bcbio_cmd= writer.bash_commands.bcbio(
        run_yaml,
        os.path.join(bcbio_dir, 'work'),
        threads=16
    )
    prep_status = executor.execute([' '.join(sample_prep)], env='local').join()
    app_logger.info('BCBio sample prep exit status: ' + str(prep_status))

    user_sample_id = clarity.get_user_sample_name(sample_id)
    if user_sample_id:
        app_logger.debug('Found user sample: ' + user_sample_id)

        with open(run_yaml, 'r') as i:
            run_config = yaml.load(i)
        run_config['fc_name'] = user_sample_id
        with open(run_yaml, 'w') as o:
            o.write(yaml.safe_dump(run_config, default_flow_style=False))

    bcbio_executor = executor.execute(
        [bcbio_cmd],
        prelim_cmds=writer.bash_commands.bcbio_env_vars(),
        job_name='bcb%s'%sample_id,
        run_id=sample_id,
        walltime=240,
        cpus=8,
        mem=64
    )
    os.chdir(original_dir)

    return bcbio_executor


def _cleanup(dataset_name):
    exit_status = 0
    #Wait for all the previous PBS step to be done writing to the folder before cleaning it up
    time.sleep(20)
    job_dir = os.path.join(cfg['jobs_dir'], dataset_name)
    cleanup_targets = [job_dir]
    intermediates_dir = cfg.get('intermediate_dir')
    if intermediates_dir:
        cleanup_targets.append(os.path.join(intermediates_dir, dataset_name))

    for t in cleanup_targets:
        app_logger.info('Cleaning up ' + t)
        try:
            shutil.rmtree(t)
        except (OSError, FileNotFoundError, NotADirectoryError) as e:
            app_logger.error(str(e))
            app_logger.warning('Could not remove: ' + t)
            exit_status += 1

    return exit_status
