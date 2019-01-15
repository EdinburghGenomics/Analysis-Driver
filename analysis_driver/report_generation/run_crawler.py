import csv, copy, json
from os.path import isfile
from collections import defaultdict, Counter
from scipy.stats.mstats import theilslopes
from egcg_core import util
from egcg_core.constants import *
from egcg_core.rest_communication import patch_entry, post_or_patch as pp
from analysis_driver.exceptions import PipelineError
from analysis_driver.reader import demultiplexing_parsers as dm, mapping_stats_parsers as mp
from analysis_driver.util.helper_functions import get_run_element_id
from .crawler import Crawler


class RunCrawler(Crawler):
    STAGE_CONVERSION = 20
    STAGE_FILTER = 30
    STAGE_MAPPING = 40

    def __init__(self, dataset, run_dir=None, stage=0):
        self.dataset = dataset
        self._populate_barcode_info_from_dataset(dataset)
        self._populate_from_lims()
        if run_dir:
            self._populate_lane_info_from_interop_metrics(run_dir)
            if stage >= self.STAGE_CONVERSION:
                self._populate_barcode_info_from_json_file(run_dir)
                self._populate_barcode_info_from_well_dup(run_dir)
                self._populate_barcode_info_from_phix_read_names(run_dir)
            if stage >= self.STAGE_FILTER:
                self._populate_barcode_info_from_fastq_filterer_files(run_dir)
                self._populate_barcode_info_from_seqtk_fqchk_files(run_dir)
            if stage >= self.STAGE_MAPPING:
                self._populate_from_mapping_stats(run_dir)
                self._populate_from_gc_bias_metrics(run_dir)

    @staticmethod
    def _update_doc_list(d, k, v):
        if k not in d:
            d[k] = []
        d[k].append(v)

    @staticmethod
    def _update_doc_set(d, k, v):
        if k not in d:
            d[k] = set()
        d[k].add(v)

    def _populate_barcode_info_from_dataset(self, dataset):
        """
        :param analysis_driver.dataset.RunDataset dataset:
        """
        self.barcodes_info = defaultdict(dict)
        self.unexpected_barcodes = {}
        self.libraries = defaultdict(dict)
        self.lanes = defaultdict(dict)
        self.run = {ELEMENT_RUN_ELEMENTS: []}
        self.projects = defaultdict(dict)

        for run_element in dataset.run_elements:
            run_element_id = '%s_%s' % (dataset.name, run_element[ELEMENT_LANE])
            barcode_info = copy.copy(run_element)
            if dataset.has_barcodes:
                run_element_id += '_' + run_element[ELEMENT_BARCODE]

            barcode_info[ELEMENT_RUN_NAME] = dataset.name
            barcode_info[ELEMENT_RUN_ELEMENT_ID] = run_element_id
            self.barcodes_info[run_element_id] = barcode_info

            # Populate the libraries
            lib = self.libraries[barcode_info[ELEMENT_LIBRARY_INTERNAL_ID]]
            lib[ELEMENT_SAMPLE_INTERNAL_ID] = barcode_info[ELEMENT_SAMPLE_INTERNAL_ID]
            lib[ELEMENT_PROJECT_ID] = barcode_info[ELEMENT_PROJECT_ID]
            lib[ELEMENT_LIBRARY_INTERNAL_ID] = barcode_info[ELEMENT_LIBRARY_INTERNAL_ID]
            self._update_doc_list(lib, k=ELEMENT_RUN_ELEMENTS, v=run_element_id)

            # Populate the projects
            proj = self.projects[barcode_info[ELEMENT_PROJECT_ID]]
            proj[ELEMENT_PROJECT_ID] = barcode_info[ELEMENT_PROJECT_ID]
            self._update_doc_set(proj, k=ELEMENT_SAMPLES, v=barcode_info[ELEMENT_SAMPLE_INTERNAL_ID])

            # Populate the lanes
            lane_id = '%s_%s' % (barcode_info[ELEMENT_RUN_NAME], run_element[ELEMENT_LANE])
            ln = self.lanes[lane_id]
            ln[ELEMENT_RUN_NAME] = barcode_info[ELEMENT_RUN_NAME]
            ln[ELEMENT_LANE_ID] = lane_id
            ln[ELEMENT_LANE_NUMBER] = run_element[ELEMENT_LANE]
            self._update_doc_list(ln, k=ELEMENT_RUN_ELEMENTS, v=run_element_id)

            # Populate the run
            self.run[ELEMENT_RUN_ELEMENTS].append(run_element_id)

            if dataset.has_barcodes:
                unknown_element_id = '%s_%s_%s' % (self.dataset.name, run_element[ELEMENT_LANE], 'unknown')
                self.barcodes_info[unknown_element_id] = {
                    ELEMENT_BARCODE: 'unknown',
                    ELEMENT_RUN_ELEMENT_ID: unknown_element_id,
                    ELEMENT_RUN_NAME: self.dataset.name,
                    ELEMENT_PROJECT_ID: 'default',
                    ELEMENT_SAMPLE_INTERNAL_ID: 'Undetermined',
                    ELEMENT_LIBRARY_INTERNAL_ID: 'Undetermined',
                    ELEMENT_LANE: run_element[ELEMENT_LANE]
                }
        for project_id in self.projects:
            self.projects[project_id][ELEMENT_SAMPLES] = list(self.projects[project_id][ELEMENT_SAMPLES])

        if dataset.has_barcodes:
            # Add the unknown to the lane
            for lane_id in self.lanes:
                lane = self.lanes[lane_id][ELEMENT_LANE_NUMBER]
                unknown = '%s_%s_%s' % (self.dataset.name, lane, 'unknown')
                self.lanes[lane_id][ELEMENT_RUN_ELEMENTS].append(unknown)

        self.run[ELEMENT_NUMBER_LANE] = len(self.lanes)

    def _populate_from_lims(self):
        for libname in self.libraries:
            self.libraries[libname].update(
                self.get_sample_information_from_lims(self.libraries[libname][ELEMENT_SAMPLE_INTERNAL_ID])
            )


    def _populate_barcode_info_from_phix_read_names(self, run_dir):
        for run_element_id, barcode_info in self.barcodes_info.items():
            if ELEMENT_BARCODE in barcode_info and barcode_info[ELEMENT_BARCODE] == 'unknown':
                read_name_file = util.find_file(
                    run_dir, 'Undetermined_S0_L00%s_phix_read_name.list' % barcode_info[ELEMENT_LANE]
                )
            else:
                read_name_file = util.find_file(
                    run_dir,
                    barcode_info[ELEMENT_PROJECT_ID],
                    barcode_info[ELEMENT_SAMPLE_INTERNAL_ID],
                    '*_S*_L00%s_phix_read_name.list' % barcode_info[ELEMENT_LANE]
                )
            if read_name_file:
                with open(read_name_file) as open_file:
                    barcode_info[ELEMENT_NB_READS_PHIX] = sum(1 for _ in open_file)
            elif barcode_info[ELEMENT_NB_READS_PASS_FILTER] == 0:
                self.info('No reads for %s, not expecting PhiX filtered file', run_element_id)
            else:
                # TODO: Not mandatory for now as there will be lots of old runs without it
                self.warning('No Phix read_name file found in %s for %s', run_dir, run_element_id)

    def _populate_barcode_info_from_seqtk_fqchk_files(self, run_dir):
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            if ELEMENT_BARCODE in barcode_info and barcode_info[ELEMENT_BARCODE] == 'unknown':
                fq_chk_files = util.find_files(
                    run_dir, 'Undetermined_S0_L00%s_R*_001.fastq.gz.fqchk' % barcode_info[ELEMENT_LANE]
                )
            else:
                fq_chk_files = util.find_files(
                    run_dir,
                    barcode_info[ELEMENT_PROJECT_ID],
                    barcode_info[ELEMENT_SAMPLE_INTERNAL_ID],
                    '*_S*_L00%s_R*_001.fastq.gz.fqchk' % barcode_info[ELEMENT_LANE]
                )

            if len(fq_chk_files) == 2:
                fq_chk_files.sort()
                nb_read, nb_base, lo_q, hi_q = dm.parse_seqtk_fqchk_file(fq_chk_files[0], q_threshold=30)
                barcode_info[ELEMENT_NB_READS_CLEANED] = nb_read
                barcode_info[ELEMENT_NB_BASE_R1_CLEANED] = nb_base
                barcode_info[ELEMENT_NB_Q30_R1_CLEANED] = hi_q
                nb_read, nb_base, lo_q, hi_q = dm.parse_seqtk_fqchk_file(fq_chk_files[1], q_threshold=30)
                barcode_info[ELEMENT_NB_BASE_R2_CLEANED] = nb_base
                barcode_info[ELEMENT_NB_Q30_R2_CLEANED] = hi_q

            elif len(fq_chk_files) == 1:
                raise PipelineError('Only one fqchk file found in %s for %s' % (run_dir, run_element_id))

            elif barcode_info[ELEMENT_NB_READS_PASS_FILTER] == 0:
                barcode_info[ELEMENT_NB_READS_CLEANED] = 0
                barcode_info[ELEMENT_NB_BASE_R1_CLEANED] = 0
                barcode_info[ELEMENT_NB_Q30_R1_CLEANED] = 0
                barcode_info[ELEMENT_NB_BASE_R2_CLEANED] = 0
                barcode_info[ELEMENT_NB_Q30_R2_CLEANED] = 0
            else:
                raise PipelineError('%s fqchk files found in %s for %s' % (len(fq_chk_files), run_dir, run_element_id))

    def _populate_barcode_info_from_fastq_filterer_files(self, run_dir):
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            if ELEMENT_BARCODE in barcode_info and barcode_info[ELEMENT_BARCODE] == 'unknown':
                fastqfilter_stats_file = util.find_file(
                    run_dir, 'Undetermined_S0_L00%s_fastqfilterer.stats' % barcode_info[ELEMENT_LANE]
                )
            else:
                fastqfilter_stats_file = util.find_file(
                    run_dir,
                    barcode_info[ELEMENT_PROJECT_ID],
                    barcode_info[ELEMENT_SAMPLE_INTERNAL_ID],
                    '*_S*_L00%s_fastqfilterer.stats' % barcode_info[ELEMENT_LANE]
                )
            if fastqfilter_stats_file:
                stats = dm.parse_fastq_filterer_stats(fastqfilter_stats_file)
                # make sure the stats can be nullable if rerun without filtering
                barcode_info[ELEMENT_TILES_FILTERED] = stats.get('remove_tiles')
                barcode_info[ELEMENT_TRIM_R1_LENGTH] = stats.get('trim_r1')
                barcode_info[ELEMENT_TRIM_R2_LENGTH] = stats.get('trim_r2')
            elif barcode_info[ELEMENT_NB_READS_PASS_FILTER] == 0:
                self.info('No reads for %s, Not expecting fastqfilter file', run_element_id)
            else:
                raise PipelineError('Cannot find fastqfilter file in %s for %s' % (run_dir, run_element_id))

    def _populate_barcode_info_from_well_dup(self, run_dir):
        welldup_files = util.find_files(run_dir, self.dataset.name + '.wellduplicate')
        if welldup_files:
            dup_per_lane = dm.parse_welldup_file(welldup_files[0])
            for run_element_id in self.barcodes_info:
                barcode_info = self.barcodes_info.get(run_element_id)
                lane = int(barcode_info.get(ELEMENT_LANE))
                if lane in dup_per_lane:
                    barcode_info[ELEMENT_LANE_PC_OPT_DUP] = dup_per_lane[lane]

    def _populate_barcode_info_from_json_file(self, run_dir):
        """Parses the Stats.json file and populates the conversion statistics and adaptor trimming details in the
        barcode_info variable. It succeeds and will replace the _populate_barcode_info_from_conversion_file() and
        _populate_barcode_info_from_adapter_file() functions."""
        json_files = util.find_files(run_dir, 'Stats', 'Stats.json')
        if json_files:
            with open(json_files[0], 'r') as json_stats:
                json_data = json.load(json_stats)

            # Call function which parses of the run elements and adapter trimmings in JSON file (previously barcodes)
            all_run_elements, unknown_run_elements, adapter_trimmed_by_id = dm.parse_json_stats(json_data, self.dataset.name)

            # To find the sum of the reads per lane, and populate barcode info
            reads_per_lane = self._populate_barcode_info_from_run_elements(all_run_elements)

            # iterating over run elements to calculate proportion of each per lane
            self._aggregate_run_element_per_lane(reads_per_lane, adapter_trimmed_by_id)

            # populating the unknown run elements array
            self._populate_unknown_elements(unknown_run_elements, reads_per_lane)
        else:
            # This file is expected. The process should stop if it is not found.
            raise FileNotFoundError()

    def _populate_barcode_info_from_run_elements(self, all_run_elements):
        # to find the sum of the reads per lane
        reads_per_lane = Counter()

        # iterate over all run elements, populating the barcode_info
        for (library, lane, barcode, clust_count,
             clust_count_pf, nb_bases, nb_bases_r1_q30, nb_bases_r2_q30) in all_run_elements:

            # incrementing the count of reads per lane
            reads_per_lane[lane] += clust_count_pf

            # retrieving barcode info, which depends on the run_element_id which differs between barcodeless and
            # multiplexed runs
            run_element_id = get_run_element_id(run_id=self.dataset.name, lane_number=lane, barcode=barcode)
            barcode_info = self.barcodes_info.get(run_element_id)

            # populating the barcode info array
            barcode_info[ELEMENT_NB_READS_SEQUENCED] = clust_count
            barcode_info[ELEMENT_NB_READS_PASS_FILTER] = clust_count_pf
            barcode_info[ELEMENT_NB_BASE_R1] = nb_bases
            barcode_info[ELEMENT_NB_BASE_R2] = nb_bases
            barcode_info[ELEMENT_NB_Q30_R1] = nb_bases_r1_q30
            barcode_info[ELEMENT_NB_Q30_R2] = nb_bases_r2_q30

        return reads_per_lane

    def _aggregate_run_element_per_lane(self, reads_per_lane, adapter_trimmed_by_id):
        # iterating over run elements to calculate proportion of each per lane
        for run_element_id in self.barcodes_info:
            barcode = self.barcodes_info[run_element_id]

            # calculating the reads per lane
            reads_for_lane = reads_per_lane.get(barcode[ELEMENT_LANE])
            if reads_for_lane > 0:
                barcode[ELEMENT_PC_READ_IN_LANE] = barcode[ELEMENT_NB_READS_PASS_FILTER] / reads_for_lane

            # populating adapter trimming data
            barcode[ELEMENT_ADAPTER_TRIM_R1] = adapter_trimmed_by_id[run_element_id]['read_1_trimmed_bases']
            barcode[ELEMENT_ADAPTER_TRIM_R2] = adapter_trimmed_by_id[run_element_id]['read_2_trimmed_bases']

    def _populate_unknown_elements(self, unknown_run_elements, reads_per_lane):
        # this helper function populates the unknown run elements array
        for lane, barcode, clust_count in unknown_run_elements:
            unknown_element_id = get_run_element_id(run_id=self.dataset.name, lane_number=lane, barcode=barcode)

            self.unexpected_barcodes[unknown_element_id] = {
                ELEMENT_RUN_ELEMENT_ID: unknown_element_id,
                ELEMENT_RUN_NAME: self.dataset.name,
                ELEMENT_LANE: lane,
                ELEMENT_PC_READ_IN_LANE: int(clust_count) / reads_per_lane.get(lane),
                ELEMENT_BARCODE: barcode,
                ELEMENT_NB_READS_PASS_FILTER: int(clust_count)
            }


    def _populate_from_mapping_stats(self, run_dir):
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            if ELEMENT_BARCODE in barcode_info and barcode_info[ELEMENT_BARCODE] == 'unknown':
                # No mapping for unassigned run element
                continue

            if barcode_info[ELEMENT_NB_READS_PASS_FILTER] == 0:
                self.info('No reads for %s, not expecting mapping stat file', run_element_id)
                continue

            fastq_file = util.find_files(
                run_dir,
                barcode_info[ELEMENT_PROJECT_ID],
                barcode_info[ELEMENT_SAMPLE_INTERNAL_ID],
                '*_S*_L00%s_R1_001.fastq.gz' % barcode_info[ELEMENT_LANE]
            )
            assert len(fastq_file) == 1
            mapping_statistics = {}
            fastq_base = fastq_file[0][:-len('_R1_001.fastq.gz')]
            samtools_stat = fastq_base + '_samtools_stats.txt'
            if isfile(samtools_stat):
                (total_reads, mapped_reads, duplicate_reads, proper_pairs) = mp.parse_samtools_stats(samtools_stat)
                mapping_statistics.update({
                    ELEMENT_NB_READS_IN_BAM: total_reads,
                    ELEMENT_NB_MAPPED_READS: mapped_reads,
                    ELEMENT_NB_DUPLICATE_READS: duplicate_reads,
                    ELEMENT_NB_PROPERLY_MAPPED: proper_pairs
                })

            samtools_depth = fastq_base + '_samtools.depth'
            if isfile(samtools_depth):
                (mean, median, sd, coverage_percentiles, bases_at_coverage,
                 genome_size, evenness) = dm.get_coverage_statistics(samtools_depth)
                coverage_statistics = {
                    ELEMENT_MEAN_COVERAGE: mean,
                    ELEMENT_MEDIAN_COVERAGE_SAMTOOLS: median,
                    ELEMENT_COVERAGE_SD: sd,
                    ELEMENT_COVERAGE_PERCENTILES: coverage_percentiles,
                    ELEMENT_BASES_AT_COVERAGE: bases_at_coverage,
                    ELEMENT_SAMPLE_GENOME_SIZE: genome_size,
                    ELEMENT_COVERAGE_EVENNESS: evenness
                }
                barcode_info[ELEMENT_COVERAGE_STATISTICS] = coverage_statistics

            picard_mark_dup_metric = fastq_base + '_markdup.metrics'
            if isfile(picard_mark_dup_metric):
                mapped_reads, dup_reads, opt_dup_reads, est_library_size = mp.parse_picard_mark_dup_metrics(picard_mark_dup_metric)
                mapping_statistics.update({
                    ELEMENT_NB_PICARD_DUP_READS: dup_reads,
                    ELEMENT_NB_PICARD_OPT_DUP_READS: opt_dup_reads,
                })
                # Estimated library size can be missing if the number of read is too low
                if est_library_size:
                    mapping_statistics[ELEMENT_PICARD_EST_LIB_SIZE] = est_library_size

            picard_insert_size_metric = fastq_base + '_insertsize.metrics'
            if isfile(picard_insert_size_metric):
                insert_types = mp.parse_picard_insert_size_metrics(picard_insert_size_metric)
                if 'FR' in insert_types:
                    mapping_statistics.update(insert_types.pop('FR'))
                if insert_types:
                    mapping_statistics[ELEMENT_NON_FR_INSERTS] = insert_types
            if mapping_statistics:
                barcode_info[ELEMENT_MAPPING_STATISTICS] = mapping_statistics

    def _populate_lane_info_from_interop_metrics(self, run_dir):
        interop_files = util.find_files(run_dir, 'interop_summary.txt')
        if interop_files:
            interop_metrics_per_lane = dm.parse_interop_summary(interop_files[0])
            for lane_id in self.lanes:
                lane_info = self.lanes.get(lane_id)
                lane = int(lane_info.get(ELEMENT_LANE_NUMBER))
                lane_info['interop_metrics'] = interop_metrics_per_lane.get(str(lane), {})

    def _populate_from_gc_bias_metrics(self, run_dir):
        for k, run_element in self.barcodes_info.items():
            if run_element.get('barcode') == 'unknown' or run_element[ELEMENT_NB_READS_PASS_FILTER] == 0:
                self.info('No reads for %s, not expecting GC bias data', run_element['run_element_id'])
                continue

            metrics_file = util.find_file(
                run_dir,
                run_element['project_id'],
                run_element['sample_id'],
                '*_S*_L00%s_gc_bias.metrics' % run_element['lane']
            )

            with open(metrics_file) as f:
                header = ''
                while not header.startswith('ACCUMULATION_LEVEL'):
                    header = f.readline()

                reader = csv.DictReader(f, header.split('\t'), delimiter='\t')
                lines = [l for l in reader]

                # gc slope
                data_points = [float(l['NORMALIZED_COVERAGE']) for l in lines if 20 <= int(l['GC']) <= 80]
                gc_slope = theilslopes(data_points)
                self.info('Calculated a GC slope of %s from %s data points', gc_slope, len(data_points))

                # deviation from normal
                total_windows = sum([int(l['WINDOWS']) for l in lines])
                # total_windows * 0.0004 gives approximately the same number of data points as 20 <= GC <= 80
                threshold = total_windows * 0.0004
                diffs = [abs(1 - float(l['NORMALIZED_COVERAGE'])) for l in lines if int(l['WINDOWS']) > threshold]
                normal_dev = sum(diffs) / len(diffs)
                self.info('Calculated a normal deviation of %s from %s data points', normal_dev, len(diffs))

                run_element['gc_bias'] = {
                    'slope': gc_slope[0],
                    'mean_deviation': normal_dev
                }

    def send_data(self):
        pp('run_elements', self.barcodes_info.values(), ELEMENT_RUN_ELEMENT_ID)
        pp('unexpected_barcodes', self.unexpected_barcodes.values(), ELEMENT_RUN_ELEMENT_ID)
        pp('lanes', self.lanes.values(), ELEMENT_LANE_ID)
        patch_entry('runs', self.run, ELEMENT_RUN_NAME, self.dataset.name)
        pp('samples', self.libraries.values(), ELEMENT_SAMPLE_INTERNAL_ID, ['run_elements'])
        pp('projects', self.projects.values(), ELEMENT_PROJECT_ID, ['samples'])
