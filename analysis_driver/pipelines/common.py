import os
import time
import shutil
from egcg_core import executor, clarity, util
from analysis_driver.util import bash_commands, bcbio_prepare_samples_cmd
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.config import output_files_config, default as cfg
from analysis_driver.report_generation.report_crawlers import SampleCrawler
from analysis_driver.transfer_data import output_sample_data, create_links_from_bcbio, prepare_sample_data
from analysis_driver.exceptions import PipelineError
from analysis_driver import quality_control as qc
from analysis_driver import segmentation
app_logger = log_cfg.get_logger('common')


class VarCallingStage(segmentation.Stage):
    @property
    def fastq_pair(self):
        return util.find_files(self.job_dir, 'merged', self.dataset.user_sample_id + '_R?.fastq.gz')

    @property
    def expected_output_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')


class MergeFastqs(VarCallingStage):
    def _run(self):
        fastq_files = prepare_sample_data(self.dataset)
        bcbio_prepare_sample(self.job_dir, self.dataset.name, fastq_files)
        return 0


class FastQC(VarCallingStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fastq_file) for fastq_file in self.fastq_pair],
            job_name='fastqc',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class BWAMem(VarCallingStage):
    def _run(self):
        return executor.execute(
            bash_commands.bwa_mem_biobambam(
                self.fastq_pair,
                self.dataset.reference_genome,
                self.expected_output_bam,
                {'ID': '1', 'SM': self.dataset.user_sample_id, 'PL': 'illumina'},
                thread=16
            ),
            job_name='bwa_mem',
            working_dir=self.job_dir,
            cpus=16,
            mem=64
        ).join()


class SamtoolsStats(VarCallingStage):
    def _run(self):
        return executor.execute(
            bash_commands.samtools_stats(
                self.expected_output_bam,
                os.path.join(self.job_dir, 'samtools_stats.txt')
            ),
            job_name='samtools',
            working_dir=self.job_dir,
            cpus=1,
            mem=8,
            log_commands=False
        ).join()


def build_bam_file_production(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    merge_fastqs = stage(MergeFastqs)
    fastqc = stage(FastQC, previous_stages=[merge_fastqs])
    bwa = stage(BWAMem, previous_stages=[MergeFastqs])
    contam = stage(qc.ContaminationCheck, previous_stages=[bwa], fastq_files=bwa.fastq_pair[:1])
    blast = stage(qc.ContaminationBlast, previous_stages=[bwa], fastq_file=bwa.fastq_pair[0])
    samtools_stat = stage(SamtoolsStats, previous_stages=[fastqc, bwa, contam, blast])
    samtools_depth = stage(qc.SamtoolsDepth, bam_file=bwa.expected_output_bam)

    return [samtools_stat, samtools_depth]


def link_results_files(sample_id, sample_dir, output_fileset):
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


def output_data(dataset, sample_dir, sample_id, dir_with_linked_files):
    exit_status = 0

    # upload the data to the rest API
    project_id = clarity.find_project_name_from_sample(sample_id)
    c = SampleCrawler(sample_id, project_id, dir_with_linked_files)
    c.send_data()

    # md5sum
    dataset.start_stage('md5sum')
    md5sum_exit_status = executor.execute(
        *[bash_commands.md5sum(os.path.join(dir_with_linked_files, f)) for f in os.listdir(dir_with_linked_files)],
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


def bcbio_prepare_sample(job_dir, sample_id, fastq_files):
    """Merge the fastq files per sample using bcbio prepare sample"""
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)
    cmd = bcbio_prepare_samples_cmd(job_dir, sample_id, fastq_files, user_sample_id=user_sample_id)

    exit_status = executor.execute(cmd, job_name='bcbio_prepare_samples', working_dir=job_dir,).join()
    sample_fastqs = util.find_files(job_dir, 'merged', user_sample_id + '_R?.fastq.gz')

    app_logger.info('bcbio_prepare_samples finished with exit status ' + str(exit_status))

    return sample_fastqs


def get_genome_version(dataset_name, species):
    genome_version = clarity.get_sample(dataset_name).udf.get('Genome Version')
    if genome_version is None:
        genome_version = cfg.query('species', species, 'default')
    reference = cfg.query('genomes', genome_version, 'fasta')
    if not reference:
        raise PipelineError('Could not find reference for species %s in sample %s ' % (species, dataset_name))
    return genome_version, reference


def get_dbsnp(genome_version):
    return cfg.query('genomes', genome_version, 'dbsnp')


def get_known_indels(genome_version):
    return cfg.query('genomes', genome_version, 'known_indels')


def cleanup(dataset_name):
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
