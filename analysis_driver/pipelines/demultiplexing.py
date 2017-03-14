import time
import shutil
from os import mkdir
from os.path import join, exists, isdir
from egcg_core import executor, clarity, util
from analysis_driver.pipelines.common import cleanup
from analysis_driver import segmentation
from analysis_driver.util import bash_commands
from analysis_driver.exceptions import SequencingRunError
from analysis_driver.quality_control.lane_duplicates import WellDuplicates
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.report_generation.report_crawlers import RunCrawler
from analysis_driver.transfer_data import output_run_data


class DemultiplexingStage(segmentation.Stage):
    @property
    def fastq_dir(self):
        return join(self.jobs_dir, 'fastq')


class Setup(DemultiplexingStage):
    def _run(self):
        self.info('Job dir: ' + self.job_dir)
        self.info('Input BCL folder: ' + self.input_dir)
        self.info('Fastq dir: ' + self.fastq_dir)

        if not isdir(self.fastq_dir):
            mkdir(self.fastq_dir)

        # Send the information about the run to the rest API
        crawler = RunCrawler(self.dataset.name, self.dataset.sample_sheet)
        crawler.send_data()

        # Need to sleep to make sure the LIMS has had enough time to update itself
        time.sleep(900)
        run_status = clarity.get_run(self.dataset.name).udf.get('Run Status')
        # TODO: catch bcl2fastq error logs instead
        if run_status != 'RunCompleted':
            self.error('Run status is \'%s\'. Stopping.', run_status)
            raise SequencingRunError(run_status)


class Bcl2FastqAndFilter(DemultiplexingStage):
    previous_stages = [Setup]

    def _run(self):
        self.info('bcl2fastq mask: ' + self.dataset.mask)  # e.g: mask = 'y150n,i6,y150n'
        bcl2fastq_exit_status = executor.execute(
            bash_commands.bcl2fastq(
                self.input_dir, self.fastq_dir, self.dataset.sample_sheet.filename, self.dataset.mask
            ),
            job_name='bcl2fastq',
            working_dir=self.job_dir,
            cpus=8,
            mem=32
        ).join()
        if bcl2fastq_exit_status:
            return bcl2fastq_exit_status

        # Copy the Samplesheet Runinfo.xml run_parameters.xml to the fastq dir
        for f in ('SampleSheet.csv', 'SampleSheet_analysis_driver.csv', 'runParameters.xml',
                  'RunInfo.xml', 'RTAConfiguration.xml'):
            shutil.copy2(join(self.input_dir, f), join(self.fastq_dir, f))
        if not exists(join(self.fastq_dir, 'InterOp')):
            shutil.copytree(join(self.input_dir, 'InterOp'), join(self.fastq_dir, 'InterOp'))

        return executor.execute(
            *[bash_commands.fastq_filterer_and_pigz_in_place(fqs) for fqs in util.find_all_fastq_pairs(self.fastq_dir)],
            job_name='fastq_filterer',
            working_dir=self.job_dir,
            cpus=18,
            mem=10
        ).join()


class _WellDuplicates(DemultiplexingStage):
    previous_stages = [Setup]

    def _run(self):
        # well duplicates
        # This could be executed at the same time as bcl2fastq but I need the fastq directory to exist
        well_dup_exec = WellDuplicates(self.dataset, self.job_dir, self.fastq_dir, self.input_dir)
        well_dup_exec.start()
        return well_dup_exec.join()


class IntegrityCheck(DemultiplexingStage):
    previous_stages = [Bcl2FastqAndFilter]

    def _run(self):
        return executor.execute(
            *[bash_commands.gzip_test(fq) for fq in util.find_all_fastqs(self.fastq_dir)],
            job_name='integrity_check',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class FastQC(DemultiplexingStage):
    previous_stages = [Bcl2FastqAndFilter]

    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fq) for fq in util.find_all_fastqs(self.fastq_dir)],
            job_name='fastqc',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class SeqtkFQChk(DemultiplexingStage):
    previous_stages = [Bcl2FastqAndFilter]

    def _run(self):
        return executor.execute(
            *[bash_commands.seqtk_fqchk(fq) for fq in util.find_all_fastqs(self.fastq_dir)],
            job_name='fqchk',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


class MD5Sum(DemultiplexingStage):
    previous_stages = [Bcl2FastqAndFilter]

    def _run(self):
        return executor.execute(
            *[bash_commands.md5sum(fq) for fq in util.find_all_fastqs(self.fastq_dir)],
            job_name='md5sum',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


class QCOutput(DemultiplexingStage):
    previous_stages = [_WellDuplicates, IntegrityCheck, FastQC, SeqtkFQChk, MD5Sum]

    def _run(self):
        # Find conversion xml file and adapter file, and send the results to the rest API
        conversion_xml = join(self.fastq_dir, 'Stats', 'ConversionStats.xml')
        adapter_trim_file = join(self.fastq_dir, 'Stats', 'AdapterTrimming.txt')

        if exists(conversion_xml) and exists(adapter_trim_file):
            self.info('Found ConversionStats and AdaptorTrimming. Sending data.')
            crawler = RunCrawler(
                self.dataset.name, self.sample_sheet, adapter_trim_file=adapter_trim_file,
                conversion_xml_file=conversion_xml, run_dir=self.fastq_dir
            )
            # TODO: review whether we need this
            json_file = join(self.fastq_dir, 'demultiplexing_results.json')
            crawler.write_json(json_file)
            crawler.send_data()
        else:
            self.error('ConversionStats or AdaptorTrimming not found.')
            return 1


class DataOutput(DemultiplexingStage):
    previous_stages = [QCOutput]

    def _run(self):
        write_versions_to_yaml(join(self.fastq_dir, 'program_versions.yaml'))
        return output_run_data(self.fastq_dir, self.dataset.name)


class Cleanup(DemultiplexingStage):
    previous_stages = [DataOutput]

    def _run(self):
        return cleanup(self.dataset.name)
