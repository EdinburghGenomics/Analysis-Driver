import os
import time
import shutil
import yaml
from analysis_driver import reader, util, executor, quality_control as qc
from analysis_driver.external_data import clarity
from analysis_driver.dataset_scanner import RunDataset, SampleDataset
from analysis_driver.exceptions import PipelineError, SequencingRunError
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver.config import output_files_config, default as cfg
from analysis_driver.quality_control.lane_duplicates import WellDuplicates
from analysis_driver.report_generation.report_crawlers import RunCrawler, SampleCrawler
from analysis_driver.transfer_data import prepare_run_data, prepare_sample_data, output_sample_data,\
    output_run_data, create_links_from_bcbio

app_logger = log_cfg.get_logger('driver')


def pipeline(dataset):

    if isinstance(dataset, RunDataset):
        return demultiplexing_pipeline(dataset)
    elif isinstance(dataset, SampleDataset):
        species = clarity.get_species_from_sample(dataset.name)
        if species is None:
            raise PipelineError('No species information found in the LIMS for ' + dataset.name)
        elif species == 'Homo sapiens':
            return variant_calling_pipeline(dataset)
        else:
            return qc_pipeline(dataset, species)
    else:
        raise AssertionError('Unexpected dataset type: ' + str(dataset))


def demultiplexing_pipeline(dataset):
    """
    :param RunDataset dataset: Dataset object
    :return: Exit status
    :rtype: int
    """
    exit_status = 0

    dataset.start_stage('transfer')
    input_run_folder = prepare_run_data(dataset)
    dataset.end_stage('transfer')

    run_id = os.path.basename(input_run_folder)
    job_dir = os.path.join(cfg['jobs_dir'], run_id)
    fastq_dir = os.path.join(job_dir, 'fastq')
    app_logger.info('Input run folder (bcl data source): ' + input_run_folder)
    app_logger.info('Fastq dir: ' + fastq_dir)
    app_logger.info('Job dir: ' + job_dir)

    run_info = reader.RunInfo(input_run_folder)
    dataset.start_stage('setup')
    reader.transform_sample_sheet(input_run_folder, remove_barcode=not run_info.mask.has_barcodes)
    sample_sheet = reader.SampleSheet(os.path.join(input_run_folder, 'SampleSheet_analysis_driver.csv'), has_barcode=run_info.mask.has_barcodes)
    if not sample_sheet.validate(run_info.mask):
        raise PipelineError('Validation failed. Check SampleSheet.csv and RunInfo.xml.')

    # Send the information about the run to the rest API
    crawler = RunCrawler(run_id, sample_sheet)
    crawler.send_data()
    dataset.end_stage('setup')

    # Need to sleep to make sure the LIMS has had enough time to update itself
    time.sleep(900)
    run_status = clarity.get_run(run_id).udf.get('Run Status')
    # TODO: catch bcl2fastq error logs instead
    if run_status != 'RunCompleted':
        app_logger.error('Run status is \'%s\'. Stopping.', run_status)
        raise SequencingRunError(run_status)

    # bcl2fastq
    mask = sample_sheet.generate_mask(run_info.mask)
    app_logger.info('bcl2fastq mask: ' + mask)  # e.g: mask = 'y150n,i6,y150n'

    dataset.start_stage('bcl2fastq')
    exit_status += executor.execute(
        util.bash_commands.bcl2fastq(input_run_folder, fastq_dir, sample_sheet.filename, mask),
        job_name='bcl2fastq',
        working_dir=job_dir,
        cpus=8,
        mem=32
    ).join()
    dataset.end_stage('bcl2fastq', exit_status)
    if exit_status:
        return exit_status

    # well duplicates
    # This could be executed at the same time as bcl2fastq but I need the fastq directory to exist
    well_dup_exec = WellDuplicates(dataset, job_dir, fastq_dir, input_run_folder)
    well_dup_exec.start()

    # Filter the adapter dimer from fastq with sickle
    dataset.start_stage('sickle_filter')
    exit_status = executor.execute(
        *[util.bash_commands.sickle_paired_end_in_place(fqs) for fqs in util.find_all_fastq_pairs(fastq_dir)],
        job_name='sickle_filter',
        working_dir=job_dir,
        cpus=1,
        mem=2
    ).join()
    dataset.end_stage('sickle_filter', exit_status)
    if exit_status:
        return exit_status


    exit_status += well_dup_exec.join()
    if exit_status:
        return exit_status

    # fastqc
    dataset.start_stage('fastqc')
    fastqc_executor = executor.execute(
        *[util.bash_commands.fastqc(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='fastqc',
        working_dir=job_dir,
        cpus=1,
        mem=2
    )

    # seqtk fqchk
    dataset.start_stage('seqtk_fqchk')
    seqtk_fqchk_executor = executor.execute(
        *[util.bash_commands.seqtk_fqchk(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='fqchk',
        working_dir=job_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    # md5sum
    dataset.start_stage('md5sum')
    md5sum_executor = executor.execute(
        *[util.bash_commands.md5sum(fq) for fq in util.find_all_fastqs(fastq_dir)],
        job_name='md5sum',
        working_dir=job_dir,
        cpus=1,
        mem=2,
        log_commands=False
    )

    fastqc_exit_status = fastqc_executor.join()
    dataset.end_stage('fastqc', fastqc_exit_status)
    seqtk_fqchk_exit_status = seqtk_fqchk_executor.join()
    dataset.end_stage('seqtk_fqchk', seqtk_fqchk_exit_status)
    md5_exit_status = md5sum_executor.join()
    dataset.end_stage('md5sum', md5_exit_status)

    # Copy the Samplesheet Runinfo.xml run_parameters.xml to the fastq dir
    for f in ['SampleSheet.csv', 'SampleSheet_analysis_driver.csv', 'runParameters.xml',
              'RunInfo.xml', 'RTAConfiguration.xml']:
        shutil.copy2(os.path.join(input_run_folder, f), os.path.join(fastq_dir, f))
    if not os.path.exists(os.path.join(fastq_dir, 'InterOp')):
        shutil.copytree(os.path.join(input_run_folder, 'InterOp'), os.path.join(fastq_dir, 'InterOp'))

    # Find conversion xml file and send the results to the rest API
    conversion_xml = os.path.join(fastq_dir, 'Stats', 'ConversionStats.xml')
    if os.path.exists(conversion_xml):
        app_logger.info('Found ConversionStats. Sending data.')
        crawler = RunCrawler(run_id, sample_sheet, conversion_xml, job_dir)
        # TODO: review whether we need this
        json_file = os.path.join(fastq_dir, 'demultiplexing_results.json')
        crawler.write_json(json_file)
        crawler.send_data()
    else:
        app_logger.error('File not found: %s', conversion_xml)
        exit_status += 1

    dataset.start_stage('data_transfer')
    transfer_exit_status = output_run_data(fastq_dir, run_id)
    dataset.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status + fastqc_exit_status + md5_exit_status

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += _cleanup(run_id)
        dataset.end_stage('cleanup', exit_status)
    return exit_status


def variant_calling_pipeline(dataset):
    """
    :param SampleDataset dataset:
    :return: Exit status
    :rtype: int
    """
    exit_status = 0
    fastq_files = prepare_sample_data(dataset)

    sample_id = dataset.name
    sample_dir = os.path.join(cfg['jobs_dir'], sample_id)
    app_logger.info('Job dir: ' + sample_dir)

    # merge fastq files
    dataset.start_stage('merge fastqs')
    fastq_pair = _bcbio_prepare_sample(sample_dir, sample_id, fastq_files)
    app_logger.debug('sample fastq files: ' + str(fastq_pair))
    dataset.end_stage('merge fastqs')

    # fastqc2
    dataset.start_stage('sample_fastqc')
    fastqc2_executor = executor.execute(
        *[util.bash_commands.fastqc(fastq_file) for fastq_file in fastq_pair],
        job_name='fastqc2',
        working_dir=sample_dir,
        cpus=1,
        mem=2
    )

    # genotype validation
    dataset.start_stage('genotype validation')
    genotype_validation = qc.GenotypeValidation(dataset, sample_dir, fastq_pair)
    genotype_validation.start()

    # species contamination check
    dataset.start_stage('species contamination check')
    species_contamination_check = qc.ContaminationCheck(dataset, sample_dir, [fastq_pair[0]])
    species_contamination_check.start()

    # bcbio
    dataset.start_stage('bcbio')
    bcbio_executor = _run_bcbio(sample_id, sample_dir, fastq_pair)

    # wait for genotype_validation, fastqc, species_contamination and bcbio to finish
    geno_valid_vcf_file, geno_valid_results = genotype_validation.join()
    app_logger.info('Written files: ' + str(geno_valid_vcf_file) + ' ' + str(geno_valid_results))
    dataset.end_stage('genotype validation', genotype_validation.exit_status)

    fastqc2_exit_status = fastqc2_executor.join()
    dataset.end_stage('sample_fastqc', fastqc2_exit_status)

    species_contamination_check.join()
    dataset.end_stage('species contamination check', species_contamination_check.exit_status)

    bcbio_exit_status = bcbio_executor.join()
    dataset.end_stage('bcbio', bcbio_exit_status)

    # sort out exit statuses
    if bcbio_exit_status:
        return bcbio_exit_status

    exit_status += fastqc2_exit_status + bcbio_exit_status

    # link the bcbio file into the final directory
    dir_with_linked_files = _link_results_files(sample_id, sample_dir, 'bcbio')

    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
    # gender detection
    vcf_file = os.path.join(dir_with_linked_files, user_sample_id + '.vcf.gz')
    gender_validation = qc.GenderValidation(dataset, sample_dir, vcf_file)
    gender_validation.start()

    # sample contamination check
    bam_file = os.path.join(dir_with_linked_files, user_sample_id + '.bam')
    sample_contam = qc.VerifyBamId(dataset, sample_dir, bam_file)
    sample_contam.start()

    #coverage statistics
    dataset.start_stage('coverage statistics')
    coverage_statistics_histogram = qc.SamtoolsDepth(dataset, sample_dir, bam_file)
    coverage_statistics_histogram.start()

    exit_status += sample_contam.join()
    exit_status += gender_validation.join()
    coverage_statistics_histogram.join()
    exit_status += coverage_statistics_histogram.exit_status
    dataset.end_stage('coverage statistics', coverage_statistics_histogram.exit_status)

    exit_status += _output_data(dataset, sample_dir, sample_id, dir_with_linked_files)

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += _cleanup(sample_id)
        dataset.end_stage('cleanup', exit_status)

    return exit_status


def qc_pipeline(dataset, species):
    exit_status = 0

    fastq_files = prepare_sample_data(dataset)

    sample_id = dataset.name
    sample_dir = os.path.join(cfg['jobs_dir'], sample_id)
    app_logger.info('Job dir: ' + sample_dir)

    reference = cfg.query('references', species.replace(' ', '_'), 'fasta')
    if not reference:
        raise PipelineError('Could not find reference for species %s in sample %s ' % (species, sample_id))

    # merge fastq files
    dataset.start_stage('merge fastqs')
    fastq_pair = _bcbio_prepare_sample(sample_dir, sample_id, fastq_files)
    app_logger.debug('sample fastq files: ' + str(fastq_pair))
    dataset.end_stage('merge fastqs')

    # fastqc2
    dataset.start_stage('sample_fastqc')
    fastqc_executor = executor.execute(
        *[util.bash_commands.fastqc(fastq_file) for fastq_file in fastq_pair],
        job_name='fastqc2',
        working_dir=sample_dir,
        cpus=1,
        mem=2
    )

    # bwa mem
    expected_output_bam = os.path.join(sample_dir, sample_id + '.bam')
    app_logger.info('align %s to %s genome found at %s', sample_id, species, reference)
    dataset.start_stage('sample_bwa')
    bwa_mem_executor = executor.execute(
        util.bash_commands.bwa_mem_samblaster(fastq_pair, reference, expected_output_bam, thread=16),
        job_name='bwa_mem',
        working_dir=sample_dir,
        cpus=12,
        mem=64
    )

    fastqc_exit_status = fastqc_executor.join()
    dataset.end_stage('sample_fastqc', fastqc_exit_status)
    bwa_exit_status = bwa_mem_executor.join()
    dataset.end_stage('sample_bwa', bwa_exit_status)

    # species contamination check
    dataset.start_stage('species contamination check')
    species_contamination_check = qc.ContaminationCheck(dataset, sample_dir, [fastq_pair[0]])
    species_contamination_check.start()
    species_contamination_check.join()
    dataset.end_stage('species contamination check', species_contamination_check.exit_status)

    dataset.start_stage('bamtools_stat')
    bamtools_stat_file = os.path.join(sample_dir, 'bamtools_stats.txt')
    bamtools_exit_status = executor.execute(
        util.bash_commands.bamtools_stats(expected_output_bam, bamtools_stat_file),
        job_name='bamtools',
        working_dir=sample_dir,
        cpus=2,
        mem=4,
        log_commands=False
    ).join()
    dataset.end_stage('bamtools_stat', bamtools_exit_status)

    dataset.start_stage('coverage statistics')
    bam_file = os.path.join(sample_dir, expected_output_bam)
    coverage_statistics_histogram = qc.SamtoolsDepth(dataset, sample_dir, bam_file)
    coverage_statistics_histogram.start()
    coverage_statistics_histogram.join()
    dataset.end_stage('coverage statistics', coverage_statistics_histogram.exit_status)

    exit_status += fastqc_exit_status + bwa_exit_status + bamtools_exit_status

    # link the bcbio file into the final directory
    dir_with_linked_files = _link_results_files(sample_id, sample_dir, 'non_human_qc')

    exit_status += _output_data(dataset, sample_dir, sample_id, dir_with_linked_files)

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += _cleanup(sample_id)
        dataset.end_stage('cleanup', exit_status)

    return exit_status


def _link_results_files(sample_id, sample_dir, output_fileset):
    dir_with_linked_files = os.path.join(sample_dir, 'linked_output_files')
    os.makedirs(dir_with_linked_files, exist_ok=True)

    # Create the links from the bcbio output to one directory
    create_links_from_bcbio(
        sample_id,
        sample_dir,
        output_files_config.query(output_fileset),
        dir_with_linked_files
    )
    return dir_with_linked_files


def _output_data(dataset, sample_dir, sample_id, dir_with_linked_files):
    exit_status = 0

    # upload the data to the rest API
    project_id = clarity.find_project_name_from_sample(sample_id)
    c = SampleCrawler(sample_id, project_id, dir_with_linked_files)
    c.send_data()

    # md5sum
    dataset.start_stage('md5sum')
    md5sum_exit_status = executor.execute(
        *[util.bash_commands.md5sum(os.path.join(dir_with_linked_files, f)) for f in os.listdir(dir_with_linked_files)],
        job_name='md5sum',
        working_dir=sample_dir,
        cpus=1,
        mem=2,
        log_commands=False
    ).join()
    dataset.end_stage('md5sum', md5sum_exit_status)

    exit_status += md5sum_exit_status

    # transfer output data
    dataset.start_stage('data_transfer')
    transfer_exit_status = output_sample_data(sample_id, dir_with_linked_files, cfg['output_dir'])
    dataset.end_stage('data_transfer', transfer_exit_status)
    exit_status += transfer_exit_status

    return exit_status


def _bcbio_prepare_sample(job_dir, sample_id, fastq_files):
    """
    Merge the fastq files per sample using bcbio prepare sample
    """
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
    cmd = util.bcbio_prepare_samples_cmd(job_dir, sample_id, fastq_files, user_sample_id=user_sample_id)

    exit_status = executor.execute(cmd, job_name='bcbio_prepare_samples', working_dir=job_dir,).join()
    sample_fastqs = util.find_files(job_dir, 'merged', user_sample_id + '_R?.fastq.gz')

    app_logger.info('bcbio_prepare_samples finished with exit status ' + str(exit_status))

    return sample_fastqs


def _run_bcbio(sample_id, sample_dir, sample_fastqs):
    run_template = os.path.join(
        os.path.dirname(__file__),
        '..', 'etc', 'bcbio_alignment_' + cfg['genome'] + '.yaml'
    )
    if not os.path.isfile(run_template):
        raise PipelineError(
            'Could not find BCBio run template ' + run_template + '. Is the correct genome set?'
        )

    original_dir = os.getcwd()
    os.chdir(sample_dir)
    app_logger.debug(str(sample_fastqs))

    app_logger.debug('Setting up sample: ' + sample_id)

    bcbio_dir = os.path.join(sample_dir, 'samples_' + sample_id + '-merged')

    sample_prep = [
        os.path.join(cfg['tools']['bcbio'], 'bin', 'bcbio_nextgen.py'),
        '-w template',
        run_template,
        bcbio_dir,
        bcbio_dir + '.csv'
    ] + sample_fastqs

    run_yaml = os.path.join(bcbio_dir, 'config', 'samples_' + sample_id + '-merged.yaml')
    bcbio_cmd = util.bash_commands.bcbio(run_yaml, os.path.join(bcbio_dir, 'work'), threads=16)

    prep_status = executor.execute(' '.join(sample_prep), env='local').join()
    app_logger.info('BCBio sample prep exit status: ' + str(prep_status))

    user_sample_id = clarity.get_user_sample_name(sample_id)
    if user_sample_id:
        new_name = user_sample_id
        app_logger.debug('Found user sample: ' + user_sample_id)
    else:
        new_name = sample_id

    with open(run_yaml, 'r') as i:
        run_config = yaml.load(i)
    run_config['fc_name'] = new_name
    with open(run_yaml, 'w') as o:
        o.write(yaml.safe_dump(run_config, default_flow_style=False))

    bcbio_executor = executor.execute(
        bcbio_cmd,
        prelim_cmds=util.bash_commands.export_env_vars(),
        job_name='bcb%s' % sample_id,
        working_dir=sample_dir,
        cpus=12,
        mem=64
    )
    os.chdir(original_dir)

    return bcbio_executor


def _cleanup(dataset_name):
    exit_status = 0
    # wait for all the previous PBS steps to be done writing to the folder before cleaning it up
    time.sleep(120)
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
