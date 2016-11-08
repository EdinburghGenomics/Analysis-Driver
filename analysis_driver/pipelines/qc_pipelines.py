import os
from egcg_core import executor, clarity, util
from analysis_driver import quality_control as qc
from analysis_driver.pipelines.common import _bcbio_prepare_sample, _link_results_files, _output_data, _cleanup
from analysis_driver.util import bash_commands
from analysis_driver.exceptions import PipelineError
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.transfer_data import prepare_sample_data

app_logger = log_cfg.get_logger('qc_pipeline')


def _bam_file_production(dataset, species):
    exit_status = 0

    fastq_files = prepare_sample_data(dataset)

    sample_id = dataset.name
    sample_dir = os.path.join(cfg['jobs_dir'], sample_id)
    app_logger.info('Job dir: ' + sample_dir)
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)

    reference = cfg.query('references', species, 'fasta')
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
        *[bash_commands.fastqc(fastq_file) for fastq_file in fastq_pair],
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
        bash_commands.bwa_mem_samblaster(
            fastq_pair,
            reference,
            expected_output_bam,
            {'ID': '1', 'SM': user_sample_id, 'PL': 'illumina'},
            thread=16
        ),
        job_name='bwa_mem',
        working_dir=sample_dir,
        cpus=16,
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
        bash_commands.bamtools_stats(expected_output_bam, bamtools_stat_file),
        job_name='bamtools',
        working_dir=sample_dir,
        cpus=1,
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

    return exit_status


def qc_pipeline(dataset, species):
    sample_id = dataset.name
    sample_dir = os.path.join(cfg['jobs_dir'], sample_id)

    exit_status = _bam_file_production(dataset, species)

    # link the bcbio file into the final directory
    dir_with_linked_files = _link_results_files(sample_id, sample_dir, 'non_human_qc')
    write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))

    exit_status += _output_data(dataset, sample_dir, sample_id, dir_with_linked_files)

    if exit_status == 0:
        dataset.start_stage('cleanup')
        exit_status += _cleanup(sample_id)
        dataset.end_stage('cleanup', exit_status)

    return exit_status
    sample_fastqs = util.find_files(job_dir, 'merged', user_sample_id + '_R?.fastq.gz')


    app_logger.info('bcbio_prepare_samples finished with exit status ' + str(exit_status))

    return sample_fastqs
