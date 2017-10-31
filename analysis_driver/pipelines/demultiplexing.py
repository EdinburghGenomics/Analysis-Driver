import shutil
from os import mkdir
from os.path import join, exists, isdir, basename, dirname

import os
from egcg_core.config import cfg
from egcg_core import executor, util, clarity
from egcg_core.constants import ELEMENT_PROJECT_ID, ELEMENT_LANE, ELEMENT_SAMPLE_INTERNAL_ID
from egcg_core.util import find_file

from analysis_driver import segmentation
from analysis_driver.util import bash_commands, find_all_fastq_pairs_for_lane, get_trim_values_for_bad_cycles
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.exceptions import SequencingRunError, AnalysisDriverError
from analysis_driver.quality_control import well_duplicates, BCLValidator, BadTileCycleDetector
from analysis_driver.report_generation import RunCrawler
from analysis_driver.transfer_data import output_data_and_archive
from analysis_driver.tool_versioning import toolset

toolset_type = 'non_human_sample_processing'
name = 'demultiplexing'


class DemultiplexingStage(segmentation.Stage):
    @property
    def fastq_dir(self):
        return join(self.job_dir, 'fastq')


class Setup(DemultiplexingStage):
    def _run(self):
        self.info('Job dir: ' + self.job_dir)
        self.info('Input BCL folder: ' + self.input_dir)
        self.info('Fastq dir: ' + self.fastq_dir)

        if not isdir(self.fastq_dir):
            mkdir(self.fastq_dir)

        # Send the information about the run to the rest API
        crawler = RunCrawler(self.dataset)
        crawler.send_data()

        b = BCLValidator(self.job_dir, self.dataset)
        b.check_bcls()

        # make sure the run is not aborted or errored before checking the bcl files
        run_status = self.dataset.lims_run.udf.get('Run Status')
        if run_status != 'RunCompleted':
            self.error('Run status is \'%s\'. Stopping.', run_status)
            raise SequencingRunError(run_status)

        invalid_bcls = b.read_invalid_files()
        if invalid_bcls:
            raise AnalysisDriverError('Some BCL files are corrupted; check %s for details', b.validation_log)

        return 0


class Bcl2Fastq(DemultiplexingStage):
    def _run(self):
        self.info('bcl2fastq mask: ' + self.dataset.mask)  # e.g: mask = 'y150n,i6,y150n'
        bcl2fastq_exit_status = executor.execute(
            bash_commands.bcl2fastq(
                self.input_dir, self.fastq_dir, self.dataset.sample_sheet_file, self.dataset.mask
            ),
            job_name='bcl2fastq',
            working_dir=self.job_dir,
            cpus=8,
            mem=32
        ).join()
        if bcl2fastq_exit_status:
            return bcl2fastq_exit_status

        # Copy the Samplesheet Runinfo.xml run_parameters.xml to the fastq dir
        for f in ('SampleSheet_analysis_driver.csv', 'runParameters.xml',
                  'RunInfo.xml', 'RTAConfiguration.xml'):
            shutil.copy2(join(self.input_dir, f), join(self.fastq_dir, f))
        if not exists(join(self.fastq_dir, 'InterOp')):
            shutil.copytree(join(self.input_dir, 'InterOp'), join(self.fastq_dir, 'InterOp'))

        return bcl2fastq_exit_status


class FastqFilter(DemultiplexingStage):
    def _run(self):

        # Find conversion xml file and adapter file, and send the results to the rest API
        conversion_xml = join(self.fastq_dir, 'Stats', 'ConversionStats.xml')
        adapter_trim_file = join(self.fastq_dir, 'Stats', 'AdapterTrimming.txt')

        if exists(conversion_xml) and exists(adapter_trim_file):
            self.info('Found ConversionStats and AdaptorTrimming. Sending data.')
            crawler = RunCrawler(
                self.dataset, adapter_trim_file=adapter_trim_file,
                conversion_xml_file=conversion_xml
            )
            crawler.send_data()

        # Assess if the lanes need filtering
        filter_lanes = {1: False, 2: False, 3: False, 4: False, 5: False, 6: False, 7: False, 8: False}
        q30_threshold = float(cfg.query('fastq_filterer', 'q30_threshold', ret_default=74))
        self.debug('Q30 threshold: %s', q30_threshold)
        for lane_metrics in self.dataset.lane_metrics:
            if q30_threshold > float(lane_metrics['pc_q30']) > 0:
                self.warning(
                    'Will apply cycle and tile filtering to lane %s: %%Q30=%s < %s',
                    lane_metrics['lane_number'],
                    lane_metrics['pc_q30'],
                    q30_threshold
                )
                filter_lanes[int(lane_metrics['lane_number'])] = True

        try:
            detector = BadTileCycleDetector(self.dataset)
            bad_tiles = detector.detect_bad_tiles()
            bad_cycles = detector.detect_bad_cycles()
        except Exception as e:
            self.error(e)
            bad_tiles = {}
            bad_cycles = {}

        cmd_list = []
        for lane in filter_lanes:
            fq_pairs = find_all_fastq_pairs_for_lane(self.fastq_dir, lane)
            kwargs = {}
            if filter_lanes[lane]:
                trim_r1, trim_r2 = get_trim_values_for_bad_cycles(bad_cycles.get(lane), self.dataset.run_info)
                kwargs = {'tiles_to_filter': bad_tiles.get(lane), 'trim_r2': trim_r2}

            for fqs in fq_pairs:
                cmd_list.append(bash_commands.fastq_filterer(fqs, **kwargs))

        return executor.execute(
            *cmd_list,
            prelim_cmds=[bash_commands.fq_filt_prelim_cmd()],
            job_name='fastq_filterer',
            working_dir=self.job_dir,
            cpus=18,
            mem=10
        ).join()


class IntegrityCheck(DemultiplexingStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.gzip_test(fq) for fq in util.find_all_fastqs(self.fastq_dir)],
            job_name='integrity_check',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class FastQC(DemultiplexingStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fq) for fq in util.find_all_fastqs(self.fastq_dir)],
            job_name='fastqc',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class SeqtkFQChk(DemultiplexingStage):
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
    def _run(self):
        # Find conversion xml file and adapter file, and send the results to the rest API
        conversion_xml = join(self.fastq_dir, 'Stats', 'ConversionStats.xml')
        adapter_trim_file = join(self.fastq_dir, 'Stats', 'AdapterTrimming.txt')

        if exists(conversion_xml) and exists(adapter_trim_file):
            self.info('Found ConversionStats and AdaptorTrimming. Sending data.')
            crawler = RunCrawler(
                self.dataset, adapter_trim_file=adapter_trim_file,
                conversion_xml_file=conversion_xml, run_dir=self.fastq_dir
            )
            crawler.send_data()
            return 0
        else:
            self.error('ConversionStats or AdaptorTrimming not found.')
            return 1


class DataOutput(DemultiplexingStage):
    def _run(self):
        toolset.write_to_yaml(join(self.fastq_dir, 'program_versions.yaml'))
        return output_data_and_archive(self.fastq_dir, join(cfg['output_dir'], self.dataset.name))


class PostDemultiplexingStage(DemultiplexingStage):

    _fastq_files = {}

    @property
    def alignment_dir(self):
        return join(self.job_dir, 'alignment')

    def fastq_pair(self, run_element):
        re_id = (
            run_element.get(ELEMENT_PROJECT_ID),
            run_element.get(ELEMENT_SAMPLE_INTERNAL_ID),
            run_element.get(ELEMENT_LANE)
        )
        if re_id not in self._fastq_files:
            fqs = util.find_fastqs(
                self.fastq_dir,
                run_element.get(ELEMENT_PROJECT_ID),
                run_element.get(ELEMENT_SAMPLE_INTERNAL_ID),
                run_element.get(ELEMENT_LANE)
            )
            fqs.sort()
            self._fastq_files[re_id] = fqs
        return self._fastq_files[re_id]

    def fastq_base(self, run_element):
        fastq_pair = self.fastq_pair(run_element)
        return fastq_pair[0][:-len('_R1_001.fastq.gz')]

    def bam_path(self, run_element):
        return join(
            self.alignment_dir,
            run_element.get(ELEMENT_PROJECT_ID),
            run_element.get(ELEMENT_SAMPLE_INTERNAL_ID),
            basename(self.fastq_base(run_element))
        ) + '.bam'


class BwaAlignMulti(PostDemultiplexingStage):
    def _run(self):
        bwa_commands = []
        self.debug('Searching for fastqs in ' + self.fastq_dir)
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                # make sure the directory where the bam file will go exists
                os.makedirs(dirname(self.bam_path(run_element)), exist_ok=True)
                # get the reference genome
                species = clarity.get_species_from_sample(run_element.get(ELEMENT_SAMPLE_INTERNAL_ID))
                default_genome_version = cfg.query('species', species, 'default')
                reference_genome = cfg.query('genomes', default_genome_version, 'fasta')
                bwa_commands.append(
                    bash_commands.bwa_mem_biobambam(
                        self.fastq_pair(run_element),
                        reference_genome,
                        self.bam_path(run_element),
                        {'ID': '1', 'SM': run_element.get(ELEMENT_SAMPLE_INTERNAL_ID), 'PL': 'illumina'},
                        thread=6
                    )
                )
        return executor.execute(
            *bwa_commands,
            job_name='bwa_mem',
            working_dir=self.job_dir,
            cpus=4,
            mem=24,
            log_commands=False
        ).join()


class SamtoolsStatsMulti(PostDemultiplexingStage):
    def _run(self):
        samtools_stats_cmds = []
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                samtools_stats_cmds.append(bash_commands.samtools_stats(
                    self.bam_path(run_element),
                    self.fastq_base(run_element) + '_samtools_stats.txt'
                ))
        return executor.execute(
            *samtools_stats_cmds,
            job_name='samtoolsstats',
            working_dir=self.job_dir,
            cpus=1,
            mem=8,
            log_commands=False
        ).join()


class SamtoolsDepthMulti(PostDemultiplexingStage):

    def _run(self):
        samtools_depth_cmds = []
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                samtools_depth_cmds.append(bash_commands.samtools_depth_command(
                    self.job_dir,
                    self.bam_path(run_element),
                    self.fastq_base(run_element) + '_samtools.depth'
                ))
        return executor.execute(
            *samtools_depth_cmds,
            job_name='samtoolsdepth',
            working_dir=self.job_dir,
            cpus=1,
            mem=6,
            log_commands=False
        ).join()


class PicardMarkDuplicateMulti(PostDemultiplexingStage):
    def _run(self):
        mark_dup_cmds = []
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                out_md_bam = self.bam_path(run_element)[:-len('.bam')] + '_markdup.bam'
                metrics_file = self.fastq_base(run_element) + '_markdup.metrics'
                mark_dup_cmds.append(bash_commands.picard_mark_dup_command(
                    self.bam_path(run_element),
                    out_md_bam,
                    metrics_file
                ))
        return executor.execute(
            *mark_dup_cmds,
            job_name='picardMD',
            working_dir=self.job_dir,
            cpus=1,
            mem=12
        ).join()


class PicardInsertSizeMulti(PostDemultiplexingStage):
    def _run(self):
        insert_size_cmds = []
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                metrics_file = self.fastq_base(run_element) + '_insertsize.metrics'
                histogram_file = self.fastq_base(run_element) + '_insertsize.pdf'
                insert_size_cmds.append(bash_commands.picard_insert_size_command(
                    self.bam_path(run_element),
                    metrics_file,
                    histogram_file
                ))
        return executor.execute(
            *insert_size_cmds,
            job_name='picardIS',
            working_dir=self.job_dir,
            cpus=1,
            mem=12
        ).join()


# Need to have a different name of class bcause the sage name is based on it.
class QCOutput2(QCOutput):
    pass


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    setup = stage(Setup)
    bcl2fastq = stage(Bcl2Fastq, previous_stages=[setup])
    fastq_filter = stage(FastqFilter, previous_stages=[bcl2fastq])
    welldups = stage(well_duplicates.WellDuplicates, run_directory=bcl2fastq.input_dir, output_directory=bcl2fastq.fastq_dir, previous_stages=[setup])
    integrity_check = stage(IntegrityCheck, previous_stages=[fastq_filter])
    fastqc = stage(FastQC, previous_stages=[fastq_filter])
    seqtk = stage(SeqtkFQChk, previous_stages=[fastq_filter])
    md5 = stage(MD5Sum, previous_stages=[fastq_filter])
    qc_output = stage(QCOutput, previous_stages=[welldups, integrity_check, fastqc, seqtk, md5])
    align_output  = stage(BwaAlignMulti, previous_stages=[qc_output])
    stats_output = stage(SamtoolsStatsMulti, previous_stages=[align_output])
    depth_output = stage(SamtoolsDepthMulti, previous_stages=[align_output])
    md_output = stage(PicardMarkDuplicateMulti, previous_stages=[align_output])
    is_output = stage(PicardInsertSizeMulti, previous_stages=[align_output])
    qc_output2 = stage(QCOutput2, previous_stages=[stats_output, depth_output, md_output, is_output])
    data_output = stage(DataOutput, previous_stages=[qc_output2])
    _cleanup = stage(Cleanup, previous_stages=[data_output])

    return _cleanup
