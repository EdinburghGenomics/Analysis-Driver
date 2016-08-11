from luigi import Parameter
from egcg_core import executor
from analysis_driver.util import bash_commands
from analysis_driver import quality_control as qc
from analysis_driver.pipeline import Stage
from analysis_driver.transfer_data import prepare_sample_data
from analysis_driver.driver import _bcbio_prepare_sample


class MergeFastqs(Stage):
    def _run(self):
        fastq_files = prepare_sample_data(self.dataset)
        _bcbio_prepare_sample(self.job_dir, self.dataset_name, fastq_files)


class Fastqc(Stage):
    previous_stages = Parameter()
    fastqs = Parameter()

    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fastq_file) for fastq_file in self.fastqs],
            job_name='fastqc2',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class CoverageStats(Stage):
    previous_stages = Parameter()
    bam_file = Parameter()

    def _run(self):
        coverage_statistics_histogram = qc.SamtoolsDepth(self.dataset, self.job_dir, self.bam_file)
        coverage_statistics_histogram.start()
        coverage_statistics_histogram.join()
        return coverage_statistics_histogram.exit_status
