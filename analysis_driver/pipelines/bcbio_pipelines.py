import os

import yaml
from egcg_core import executor, clarity
from analysis_driver import quality_control as qc
from analysis_driver.exceptions import PipelineError
from analysis_driver.pipelines.common import bcbio_prepare_sample, link_results_files, output_data, cleanup
from analysis_driver.util import bash_commands
from analysis_driver.dataset_scanner import SampleDataset
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.transfer_data import prepare_sample_data

app_logger = log_cfg.get_logger('bcbio_pipelines')


def bcbio_var_calling_pipeline(dataset, genome_version, analysis_type):
    """
    :param SampleDataset dataset:
    :param analysis_type:
    :param genome_version:
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
    fastq_pair = bcbio_prepare_sample(sample_dir, sample_id, fastq_files)
    app_logger.debug('sample fastq files: ' + str(fastq_pair))
    dataset.end_stage('merge fastqs')

    # fastqc2
    dataset.start_stage('sample_fastqc')
    fastqc2_executor = executor.execute(
        *[bash_commands.fastqc(fastq_file) for fastq_file in fastq_pair],
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

    # blast contamination check
    blast_contamination_check = qc.ContaminationBlast(dataset, sample_dir, [fastq_pair[0]])
    blast_contamination_check.start()

    # bcbio
    dataset.start_stage('bcbio')
    bcbio_executor = _run_bcbio(sample_id, sample_dir, fastq_pair,
                                genome_version=genome_version, analysis_type=analysis_type)

    # wait for genotype_validation, fastqc, species_contamination and bcbio to finish
    geno_valid_vcf_file, geno_valid_results = genotype_validation.join()
    app_logger.info('Written files: ' + str(geno_valid_vcf_file) + ' ' + str(geno_valid_results))
    dataset.end_stage('genotype validation', genotype_validation.exit_status)

    fastqc2_exit_status = fastqc2_executor.join()
    dataset.end_stage('sample_fastqc', fastqc2_exit_status)

    blast_contamination_check.join()
    species_contamination_check.join()
    contam_check_status = species_contamination_check.exit_status + blast_contamination_check.exit_status
    dataset.end_stage('species contamination check', contam_check_status)

    bcbio_exit_status = bcbio_executor.join()
    dataset.end_stage('bcbio', bcbio_exit_status)

    # sort out exit statuses
    if bcbio_exit_status:
        return bcbio_exit_status

    exit_status += fastqc2_exit_status + bcbio_exit_status

    # link the bcbio file into the final directory
    dir_with_linked_files = link_results_files(sample_id, sample_dir, 'bcbio')

    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
    # gender detection
    vcf_file = os.path.join(dir_with_linked_files, user_sample_id + '.vcf.gz')
    gender_validation = qc.GenderValidation(dataset, sample_dir, vcf_file)
    gender_validation.start()

    # vcf stats
    vcf_stats = qc.VCFStats(dataset, sample_dir, vcf_file)
    vcf_stats.start()

    # sample contamination check
    bam_file = os.path.join(dir_with_linked_files, user_sample_id + '.bam')
    sample_contam = qc.VerifyBamId(dataset, sample_dir, bam_file)
    sample_contam.start()

    # coverage statistics
    dataset.start_stage('coverage statistics')
    coverage_statistics_histogram = qc.SamtoolsDepth(dataset, sample_dir, bam_file)
    coverage_statistics_histogram.start()

    exit_status += sample_contam.join()
    exit_status += gender_validation.join()
    exit_status += vcf_stats.join()
    coverage_statistics_histogram.join()
    exit_status += coverage_statistics_histogram.exit_status
    dataset.end_stage('coverage statistics', coverage_statistics_histogram.exit_status)

    write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))
    exit_status += output_data(dataset, sample_dir, sample_id, dir_with_linked_files)

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += cleanup(sample_id)
        dataset.end_stage('cleanup', exit_status)

    return exit_status


def _run_bcbio(sample_id, sample_dir, sample_fastqs, genome_version, analysis_type):
    if not genome_version:
        genome_version = cfg['genome']
    if not analysis_type:
        analysis_type = 'gatk'
    elif analysis_type.endswith('gatk'):
        analysis_type = 'gatk'
    elif analysis_type.endswith('freebayes'):
        analysis_type = 'freebayes'
    else:
        raise PipelineError('Unknown Analysis type %s' % analysis_type)

    run_template = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..', '..', 'etc', 'bcbio_alignment_%s_%s.yaml' % (genome_version, analysis_type)
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
    bcbio_cmd = bash_commands.bcbio(run_yaml, os.path.join(bcbio_dir, 'work'), threads=16)

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
        prelim_cmds=bash_commands.export_env_vars(),
        job_name='bcb%s' % sample_id,
        working_dir=sample_dir,
        cpus=12,
        mem=64
    )
    os.chdir(original_dir)

    return bcbio_executor
