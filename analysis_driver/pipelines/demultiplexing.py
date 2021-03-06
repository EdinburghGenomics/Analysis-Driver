import shutil
import time
from os import makedirs
from os.path import join, exists, dirname, basename
from egcg_core.config import cfg
from egcg_core import executor, util, rest_communication, constants as c
from analysis_driver import segmentation
from analysis_driver.reader.run_info import Reads
from analysis_driver.util import bash_commands, find_all_fastq_pairs_for_lane, get_trim_values_for_bad_cycles
from analysis_driver.pipelines import Pipeline, rapid
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.exceptions import SequencingRunError, AnalysisDriverError
from analysis_driver.quality_control import well_duplicates, interop_metrics, BCLValidator, BadTileCycleDetector
from analysis_driver.report_generation import RunCrawler
from analysis_driver.transfer_data import output_data_and_archive
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.helper_functions import merge_lane_directories


class DemultiplexingStage(segmentation.Stage):
    @property
    def fastq_dir(self):
        return join(self.job_dir, 'fastq')


class Setup(DemultiplexingStage):
    def _run(self):
        self.info('Input BCL folder: ' + self.input_dir)
        self.info('Fastq dir: ' + self.fastq_dir)

        makedirs(self.fastq_dir, exist_ok=True)

        # Send the information about the run to the rest API
        crawler = RunCrawler(self.dataset)
        crawler.send_data()

        b = BCLValidator(dataset=self.dataset)
        b.check_bcls()

        # Run is finished now so get the interop summary
        summary_file = join(self.fastq_dir, 'interop_summary.txt')
        cmd = '%s %s > %s' % (toolset['interop_summary'], self.input_dir, summary_file)
        summary_exit_status = executor.execute(
            cmd,
            job_name='summary',
            working_dir=self.job_dir,
            cpus=1,
            mem=4
        ).join()

        if summary_exit_status != 0:
            # Failure to generate the metrics is a sign of file corruptions which will be picked up after
            self.warning('Interop Summary generation failed')
        else:
            # Send the results of interop summary to the rest API
            RunCrawler(self.dataset, run_dir=self.fastq_dir).send_data()

        # make sure the run is not aborted or errored before checking the bcl files
        run_status = self.dataset.run_status.get('run_status')
        if run_status != 'RunCompleted':
            self.error('Run status is \'%s\'. Stopping.', run_status)
            raise SequencingRunError(run_status)

        invalid_bcls = b.read_invalid_files()
        if invalid_bcls:
            raise AnalysisDriverError('Some BCL files are corrupted; check %s for details', b.validation_log)

        return 0


class Bcl2Fastq(DemultiplexingStage):
    def _run(self):
        lanes = range(1, self.dataset.number_of_lane + 1)
        masks = [self.dataset.mask_per_lane(l) for l in lanes]
        self.info('bcl2fastq mask: ' + ', '.join(('lane %s: %s' % (l, m) for l, m in zip(lanes, masks))))

        bcl2fastq_exit_status = executor.execute(
            *(bash_commands.bcl2fastq_per_lane(self.input_dir, self.fastq_dir, self.dataset.sample_sheet_file_for_lane,
                                               masks, lanes)),
            job_name='bcl2fastq',
            working_dir=self.job_dir,
            cpus=8,
            mem=32
        ).join()
        if bcl2fastq_exit_status:
            return bcl2fastq_exit_status

        # Merge the lane directories
        merge_lane_directories(self.fastq_dir, self.dataset.run_elements)

        # Copy the Samplesheet Runinfo.xml run_parameters.xml to the fastq dir
        for f in ('runParameters.xml', 'RunInfo.xml', 'RTAConfiguration.xml'):
            shutil.copy2(join(self.input_dir, f), join(self.fastq_dir, f))
        for lane in range(1, 9):
            f = basename(self.dataset.sample_sheet_file_for_lane(lane))
            shutil.copy2(self.dataset.sample_sheet_file_for_lane(lane), join(self.fastq_dir, f))
        if not exists(join(self.fastq_dir, 'InterOp')):
            shutil.copytree(join(self.input_dir, 'InterOp'), join(self.fastq_dir, 'InterOp'))

        return bcl2fastq_exit_status


class PhixDetection(DemultiplexingStage):
    def _run(self):
        fasta = cfg.get('genomes_dir', '') + \
                rest_communication.get_document('genomes', where={'assembly_name': 'phix174'})['data_files']['fasta']
        cmds = []
        for fq1, fq2 in util.find_all_fastq_pairs(self.fastq_dir):
            read_name_list = fq1[:-len('_R1_001.fastq.gz')] + '_phix_read_name.list'
            cmds.append(bash_commands.bwa_mem_phix(fq1, read_name_list, fasta))

        exit_status = executor.execute(
            *cmds,
            job_name='phix_detection',
            working_dir=self.job_dir,
            cpus=16,
            mem=10,
            log_commands=False
        ).join()
        if exit_status == 0:
            # Send the results of BCL2fastq to the rest API
            crawler = RunCrawler(self.dataset, run_dir=self.fastq_dir, stage=RunCrawler.STAGE_CONVERSION)
            crawler.send_data()
        return exit_status


class FastqFilter(DemultiplexingStage):
    def _run(self):
        # Assess if the lanes need filtering
        q30_threshold = float(cfg.query('fastq_filterer', 'q30_threshold', ret_default=74))
        self.info('Q30 threshold: %s', q30_threshold)
        filter_lanes = {1: False, 2: False, 3: False, 4: False, 5: False, 6: False, 7: False, 8: False}
        for lane in self.dataset.lane_metrics:
            if q30_threshold > float(util.query_dict(lane, 'aggregated.pc_q30', ret_default=0)) > 0:
                self.warning(
                    'Will apply cycle and tile filtering to lane %s: %%Q30=%s < %s',
                    lane['lane_number'],
                    lane['aggregated']['pc_q30'],
                    q30_threshold
                )
                filter_lanes[int(lane['lane_number'])] = True

        try:
            detector = BadTileCycleDetector(self.dataset)
            bad_tiles = detector.detect_bad_tiles()
            bad_cycles = detector.detect_bad_cycles()
        except Exception as e:
            self.error(e)
            bad_tiles = {}
            bad_cycles = {}

        cmds = []
        for lane in filter_lanes:
            fq_pairs = find_all_fastq_pairs_for_lane(self.fastq_dir, lane)
            kwargs = {}
            if filter_lanes[lane]:
                trim_r1, trim_r2 = get_trim_values_for_bad_cycles(bad_cycles.get(lane), self.dataset.run_info)
                kwargs = {'rm_tiles': bad_tiles.get(lane), 'trim_r2': trim_r2}

            for fqs in fq_pairs:
                read_name_list = fqs[0][:-len('_R1_001.fastq.gz')] + '_phix_read_name.list'
                cmds.append(bash_commands.fastq_filterer(fqs, rm_reads=read_name_list, **kwargs))

        return executor.execute(
            *cmds,
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
    run_crawler_stage = segmentation.Parameter()

    def _run(self):
        crawler = RunCrawler(self.dataset, run_dir=self.fastq_dir, stage=self.run_crawler_stage)
        crawler.send_data()
        return 0


class DataOutput(DemultiplexingStage):
    def _run(self):
        toolset.write_to_yaml(join(self.fastq_dir, 'program_versions.yaml'))
        return output_data_and_archive(self.fastq_dir, join(cfg['run']['output_dir'], self.dataset.name))


class WaitForRead2(DemultiplexingStage):
    """Wait for the first 50 cycles of read 2 to be extracted."""

    def _run(self):
        # Full first read, the indexes, and 50 bases of read 2
        required_number_of_cycles = sum([
            Reads.num_cycles(self.dataset.run_info.reads.upstream_read),
            sum(self.dataset.run_info.reads.index_lengths),
            50
        ])

        # get the last cycle extracted
        current_cycle = interop_metrics.get_last_cycles_with_existing_bcls(self.dataset.input_dir)
        while current_cycle < required_number_of_cycles and self.dataset.is_sequencing():
            time.sleep(1200)
            current_cycle = interop_metrics.get_last_cycles_with_existing_bcls(self.dataset.input_dir)

        # make sure the run is not aborted or errored before continuing with the rest of the pipeline
        run_status = self.dataset.run_status.get('run_status')
        if run_status not in ['RunCompleted', 'RunStarted', 'RunPaused']:
            self.error('Run status is \'%s\'. Stopping.', run_status)
            raise SequencingRunError(run_status)

        return 0


class PartialRun(DemultiplexingStage):
    @property
    def fastq_intermediate_dir(self):
        return join(self.job_dir, 'fastq_intermediate')


class Bcl2FastqPartialRun(PartialRun):
    def _run(self):
        lanes = range(1, self.dataset.number_of_lane + 1)
        masks = [self.dataset.mask_per_lane(l).replace('y150n', 'y50n101') for l in lanes]
        self.info('bcl2fastq mask: ' + ', '.join(('lane %s: %s' % (l, m) for l,m in zip(lanes, masks))))

        bcl2fastq_exit_status = executor.execute(
            *(bash_commands.bcl2fastq_per_lane(self.input_dir, self.fastq_intermediate_dir,
                                               self.dataset.sample_sheet_file_for_lane, masks, lanes)),
            job_name='bcl2fastq_intermediate',
            working_dir=self.job_dir,
            cpus=8,
            mem=32
        ).join()
        if bcl2fastq_exit_status:
            return bcl2fastq_exit_status

        # Merge the lane directories
        merge_lane_directories(self.fastq_intermediate_dir, self.dataset.run_elements)
        return bcl2fastq_exit_status


class EarlyFastqFilter(PartialRun):
    def _run(self):
        cmds = []
        for lane in range(1, 9):
            fq_pairs = find_all_fastq_pairs_for_lane(self.fastq_intermediate_dir, lane)
            for fqs in fq_pairs:
                cmds.append(bash_commands.fastq_filterer(fqs))

        return executor.execute(
            *cmds,
            prelim_cmds=[bash_commands.fq_filt_prelim_cmd()],
            job_name='fastq_filterer',
            working_dir=self.job_dir,
            cpus=10,
            mem=6
        ).join()


class PostDemultiplexingStage(PartialRun):
    _fastq_files = {}

    @property
    def alignment_dir(self):
        return join(self.job_dir, 'alignment')

    def fastq_pair(self, run_element):
        re_id = (
            run_element.get(c.ELEMENT_PROJECT_ID),
            run_element.get(c.ELEMENT_SAMPLE_INTERNAL_ID),
            run_element.get(c.ELEMENT_LANE)
        )
        if re_id not in self._fastq_files:
            fqs = util.find_fastqs(
                self.fastq_intermediate_dir,
                run_element.get(c.ELEMENT_PROJECT_ID),
                run_element.get(c.ELEMENT_SAMPLE_INTERNAL_ID),
                run_element.get(c.ELEMENT_LANE)
            )
            fqs.sort()
            self._fastq_files[re_id] = fqs
        return self._fastq_files[re_id]

    def intermediate_fastq_base(self, run_element):
        fastq_pair = self.fastq_pair(run_element)
        return fastq_pair[0][:-len('_R1_001.fastq.gz')]

    def final_fastq_base(self, run_element):
        return join(
            self.fastq_dir,
            run_element.get(c.ELEMENT_PROJECT_ID),
            run_element.get(c.ELEMENT_SAMPLE_INTERNAL_ID),
            basename(self.intermediate_fastq_base(run_element))
        )

    def bam_path(self, run_element):
        return join(
            self.alignment_dir,
            run_element.get(c.ELEMENT_PROJECT_ID),
            run_element.get(c.ELEMENT_SAMPLE_INTERNAL_ID),
            basename(self.intermediate_fastq_base(run_element))
        ) + '.bam'


class BwaAlignMulti(PostDemultiplexingStage):
    def _run(self):
        bwa_commands = []
        self.debug('Searching for fastqs in ' + self.fastq_intermediate_dir)
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                # make sure the directory where the bam and qc files will go exists
                makedirs(dirname(self.bam_path(run_element)), exist_ok=True)
                makedirs(dirname(self.final_fastq_base(run_element)), exist_ok=True)
                bwa_commands.append(
                    bash_commands.bwa_mem_biobambam(
                        self.fastq_pair(run_element),
                        self.dataset.reference_genome(run_element),
                        self.bam_path(run_element),
                        {'ID': '1', 'SM': run_element.get(c.ELEMENT_SAMPLE_INTERNAL_ID), 'PL': 'illumina'},
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
                    self.final_fastq_base(run_element) + '_samtools_stats.txt'
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
                    self.final_fastq_base(run_element) + '_samtools.depth'
                ))
        return executor.execute(
            *samtools_depth_cmds,
            job_name='samtoolsdepth',
            working_dir=self.job_dir,
            cpus=1,
            mem=6,
            log_commands=False
        ).join()


class PicardGCBias(PostDemultiplexingStage):
    def _run(self):
        cmds = []
        for r in self.dataset.run_elements:
            if self.fastq_pair(r):
                metrics_basename = self.final_fastq_base(r) + '_gc_bias'
                cmds.append(
                    bash_commands.picard_gc_bias(
                        self.bam_path(r),
                        metrics_basename + '.metrics',
                        metrics_basename + '_summary.metrics',
                        metrics_basename + '.pdf',
                        self.dataset.reference_genome(r)
                    )
                )

        return executor.execute(
            *cmds,
            job_name='gc_bias',
            working_dir=self.job_dir,
            cpus=1,
            mem=4
        ).join()


class PicardMarkDuplicateMulti(PostDemultiplexingStage):
    def _run(self):
        mark_dup_cmds = []
        for run_element in self.dataset.run_elements:
            if self.fastq_pair(run_element):
                out_md_bam = self.bam_path(run_element)[:-len('.bam')] + '_markdup.bam'
                metrics_file = self.final_fastq_base(run_element) + '_markdup.metrics'
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
                metrics_file = self.final_fastq_base(run_element) + '_insertsize.metrics'
                histogram_file = self.final_fastq_base(run_element) + '_insertsize.pdf'
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


class RunReview(DemultiplexingStage):
    def _run(self):
        rest_communication.post_entry(
            'actions',
            {'action_type': 'automatic_run_review', 'run_id': self.dataset.name},
            use_data=True
        )
        return 0


class Demultiplexing(Pipeline):
    toolset_type = 'run_processing'

    def build(self):
        wait_for_read2 = self.stage(WaitForRead2)
        bcl2fastq_half_run = self.stage(Bcl2FastqPartialRun, previous_stages=[wait_for_read2])
        early_fastq_filter = self.stage(EarlyFastqFilter, previous_stages=[bcl2fastq_half_run])
        align_output = self.stage(BwaAlignMulti, previous_stages=[early_fastq_filter])
        stats_output = self.stage(SamtoolsStatsMulti, previous_stages=[align_output])
        depth_output = self.stage(SamtoolsDepthMulti, previous_stages=[align_output])
        md_output = self.stage(PicardMarkDuplicateMulti, previous_stages=[align_output])
        is_output = self.stage(PicardInsertSizeMulti, previous_stages=[align_output])
        gc_bias = self.stage(PicardGCBias, previous_stages=[align_output])

        setup = self.stage(Setup)
        bcl2fastq = self.stage(Bcl2Fastq, previous_stages=[setup])
        phix_detection = self.stage(PhixDetection, previous_stages=[bcl2fastq])
        fastq_filter = self.stage(FastqFilter, previous_stages=[phix_detection])
        welldups = self.stage(well_duplicates.WellDuplicates, run_directory=bcl2fastq.input_dir,
                              output_directory=bcl2fastq.fastq_dir, previous_stages=[setup])
        integrity_check = self.stage(IntegrityCheck, previous_stages=[fastq_filter])
        fastqc = self.stage(FastQC, previous_stages=[fastq_filter])
        seqtk = self.stage(SeqtkFQChk, previous_stages=[fastq_filter])
        md5 = self.stage(MD5Sum, previous_stages=[fastq_filter])
        qc_output = self.stage(QCOutput, stage_name='qcoutput1', run_crawler_stage=RunCrawler.STAGE_FILTER,
                               previous_stages=[welldups, integrity_check, fastqc, seqtk, md5])

        qc_output2 = self.stage(QCOutput, stage_name='qcoutput2', run_crawler_stage=RunCrawler.STAGE_MAPPING,
                                previous_stages=[qc_output, stats_output, depth_output, md_output, is_output, gc_bias])
        data_output = self.stage(DataOutput, previous_stages=[qc_output2])
        final_checkpoint = [data_output]
        if self.dataset.rapid_samples_by_lane:
            final_checkpoint.append(rapid.Rapid(self.dataset).build(setup))

        cleanup = self.stage(Cleanup, previous_stages=final_checkpoint)
        review = self.stage(RunReview, previous_stages=[cleanup])
        return review
