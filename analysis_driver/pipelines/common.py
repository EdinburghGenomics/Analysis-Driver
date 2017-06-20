import os
import csv
import time
import shutil
from luigi import Parameter
from egcg_core import executor, clarity, util
from egcg_core.constants import ELEMENT_PROJECT_ID, ELEMENT_LANE, ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_NAME
from analysis_driver import segmentation, quality_control as qc
from analysis_driver.util import bash_commands
from analysis_driver.config import default as cfg, OutputFileConfiguration
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.report_generation import SampleCrawler
from analysis_driver.transfer_data import output_data_and_archive, create_output_links
from analysis_driver.exceptions import PipelineError


class VarCallingStage(segmentation.Stage):
    @property
    def fq_pattern(self):
        return os.path.join(self.job_dir, 'merged', self.dataset.user_sample_id + '_R?.fastq.gz')

    @property
    def exp_bam_path(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')


class FastQC(VarCallingStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(f) for f in util.find_files(self.fq_pattern)],
            job_name='fastqc',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class BWAMem(VarCallingStage):
    def _run(self):
        return executor.execute(
            bash_commands.bwa_mem_biobambam(
                util.find_files(self.fq_pattern),
                self.dataset.reference_genome,
                self.exp_bam_path,
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
                self.exp_bam_path,
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
    bwa = stage(BWAMem, previous_stages=[merge_fastqs])
    contam = stage(qc.ContaminationCheck, previous_stages=[bwa], fq_pattern=bwa.fq_pattern)
    blast = stage(qc.ContaminationBlast, previous_stages=[bwa], fastq_file=bwa.fq_pattern.replace('?', '1'))
    samtools_stat = stage(SamtoolsStats, previous_stages=[fastqc, bwa, contam, blast])
    samtools_depth = stage(qc.SamtoolsDepth, bam_file=bwa.exp_bam_path, previous_stages=[bwa])

    return [samtools_stat, samtools_depth]


class SampleDataOutput(segmentation.Stage):
    output_fileset = Parameter()

    @property
    def output_cfg(self):
        return OutputFileConfiguration(self.output_fileset)

    def _run(self):
        dir_with_linked_files = self.link_files()
        write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))
        return self.output_data(dir_with_linked_files)

    def link_files(self):
        dir_with_linked_files = os.path.join(self.job_dir, 'linked_output_files')
        os.makedirs(dir_with_linked_files, exist_ok=True)

        # Create the links from the bcbio output to one directory
        create_output_links(self.dataset.name, self.job_dir, self.output_cfg, dir_with_linked_files)
        return dir_with_linked_files

    def output_data(self, dir_with_linked_files):
        # upload the data to the rest API
        project_id = clarity.find_project_name_from_sample(self.dataset.name)
        c = SampleCrawler(self.dataset.name, project_id, self.job_dir, self.output_cfg)
        c.send_data()

        # md5sum
        self.dataset.start_stage('md5sum')
        md5sum_exit_status = executor.execute(
            *[bash_commands.md5sum(os.path.join(dir_with_linked_files, f)) for f in os.listdir(dir_with_linked_files)],
            job_name='md5sum',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()
        self.dataset.end_stage('md5sum', md5sum_exit_status)

        # transfer output data
        output_dir = os.path.join(cfg['output_dir'], project_id, self.dataset.name)
        transfer_exit_status = output_data_and_archive(
            dir_with_linked_files.rstrip('/') + '/',
            output_dir.rstrip('/')
        )
        return md5sum_exit_status + transfer_exit_status


class MergeFastqs(VarCallingStage):
    def _find_fastqs_for_run_element(self, run_element):
        local_fastq_dir = os.path.join(cfg['input_dir'], run_element.get(ELEMENT_RUN_NAME))
        self.debug('Searching for fastqs in ' + local_fastq_dir)
        return util.find_fastqs(
            local_fastq_dir,
            run_element.get(ELEMENT_PROJECT_ID),
            self.dataset.name,
            run_element.get(ELEMENT_LANE)
        )

    def find_fastqs_for_sample(self):
        self.debug('Preparing dataset %s (%s)', self.dataset.name, self.dataset.dataset_status)
        fastqs = []
        for run_element in self.dataset.run_elements:
            if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                fastqs.extend(self._find_fastqs_for_sample(run_element))
        return fastqs

    def _write_bcbio_csv(self, fastqs):
        """Write out a simple csv mapping fastq files to a sample id."""
        csv_file = os.path.join(self.job_dir, 'samples_' + self.dataset.name + '.csv')
        self.info('Writing BCBio sample csv ' + csv_file)

        with open(csv_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['samplename', 'description'])
            for fq in fastqs:
                writer.writerow([fq, self.dataset.user_sample_id])

        return csv_file

    def _run(self):
        """Merge the fastq files per sample using bcbio prepare sample"""
        fastq_files = self.prepare_sample_data()
        bcbio_csv_file = self._write_bcbio_csv(fastq_files)
        self.info('Setting up BCBio samples from ' + bcbio_csv_file)
        cmd = bash_commands.bcbio_prepare_samples(
            self.job_dir, bcbio_csv_file
        )

        exit_status = executor.execute(cmd, job_name='bcbio_prepare_samples', working_dir=self.job_dir).join()
        sample_fastqs = util.find_files(self.job_dir, 'merged', self.dataset.user_sample_id + '_R?.fastq.gz')

        self.info(
            'bcbio_prepare_samples finished with exit status %s. Merged fastqs: %s',
            exit_status,
            sample_fastqs
        )
        return 0


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


class Cleanup(segmentation.Stage):
    def _run(self):
        # wait for all the previous PBS steps to be done writing to the folder before cleaning it up
        time.sleep(120)
        job_dir = os.path.join(cfg['jobs_dir'], self.dataset.name)

        self.info('Cleaning up job dir %s', job_dir)
        try:
            shutil.rmtree(job_dir)
            return 0
        except (OSError, FileNotFoundError, NotADirectoryError) as e:
            self.error('Could not remove job dir: %s', e)
            return 1
