import copy
from collections import defaultdict, Counter
from egcg_core import util
from egcg_core.constants import *
from egcg_core.rest_communication import post_or_patch as pp
from analysis_driver.reader import demultiplexing_parsers as dm
from analysis_driver.exceptions import PipelineError
from .crawler import Crawler


class RunCrawler(Crawler):
    def __init__(self, dataset, adapter_trim_file=None, conversion_xml_file=None, run_dir=None):
        self.dataset = dataset
        self.adapter_trim_file = adapter_trim_file
        self._populate_barcode_info_from_dataset(dataset)
        self._populate_from_lims()
        if adapter_trim_file:
            self._populate_barcode_info_from_adapter_file(adapter_trim_file)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)
        if run_dir:
            self._populate_barcode_info_from_seqtk_fqchk_files(run_dir)
            welldup_files = util.find_files(run_dir, self.dataset.name + '.wellduplicate')
            if welldup_files:
                self._populate_barcode_info_from_well_dup(welldup_files[0])

            self._populate_barcode_info_from_fastq_filterer_files(run_dir)

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
        self.run = {ELEMENT_RUN_NAME: self.dataset.name, ELEMENT_RUN_ELEMENTS: []}
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

    def _run_sample_lane_to_barcode(self, adapters_trimmed_by_id):
        run_element_adapters_trimmed = {}
        for adapter_id in adapters_trimmed_by_id:
            run_element_id = None
            run_id, sample_id, lane = adapter_id
            if self.dataset.has_barcodes:
                for i in self.barcodes_info:
                    if self.barcodes_info[i][ELEMENT_RUN_NAME] == run_id \
                            and self.barcodes_info[i][ELEMENT_LANE] == lane:
                        if self.barcodes_info[i][ELEMENT_SAMPLE_INTERNAL_ID] == sample_id:
                            run_element_id = self.barcodes_info[i][ELEMENT_RUN_ELEMENT_ID]
                        elif self.barcodes_info[i][ELEMENT_SAMPLE_INTERNAL_ID] == 'Undetermined' and sample_id == 'unknown':
                            run_element_id = self.barcodes_info[i][ELEMENT_RUN_ELEMENT_ID]
            else:
                run_element_id = '%s_%s' % (run_id, lane)
            run_element_adapters_trimmed[run_element_id] = adapters_trimmed_by_id[adapter_id]
        return run_element_adapters_trimmed

    def _populate_barcode_info_from_adapter_file(self, adapter_trim_file):
        parsed_trimmed_adapters = dm.parse_adapter_trim_file(adapter_trim_file, self.dataset.name)
        run_element_adapters_trimmed = self._run_sample_lane_to_barcode(parsed_trimmed_adapters)
        for run_element_id in self.barcodes_info:
            self.barcodes_info[run_element_id][ELEMENT_ADAPTER_TRIM_R1] = run_element_adapters_trimmed[run_element_id]['read_1_trimmed_bases']
            self.barcodes_info[run_element_id][ELEMENT_ADAPTER_TRIM_R2] = run_element_adapters_trimmed[run_element_id]['read_2_trimmed_bases']

    def _populate_barcode_info_from_seqtk_fqchk_files(self, run_dir):
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            if ELEMENT_BARCODE in barcode_info and barcode_info[ELEMENT_BARCODE] == 'unknown':
                fq_chk_files = util.find_files(
                    run_dir,
                    'Undetermined_S0_L00%s_R*_001.fastq.gz.fqchk' % barcode_info[ELEMENT_LANE]
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
                    run_dir,
                    'Undetermined_S0_L00%s_fastqfilterer.stats' % barcode_info[ELEMENT_LANE]
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

    def _populate_barcode_info_from_well_dup(self, welldup_file):
        dup_per_lane = dm.parse_welldup_file(welldup_file)
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            lane = int(barcode_info.get(ELEMENT_LANE))
            if lane in dup_per_lane:
                barcode_info[ELEMENT_LANE_PC_OPT_DUP] = dup_per_lane[lane]

    def _populate_barcode_info_from_conversion_file(self, conversion_xml):
        all_barcodes, top_unknown_barcodes, all_barcodeless = dm.parse_conversion_stats(conversion_xml, self.dataset.has_barcodes)
        reads_per_lane = Counter()
        if self.dataset.has_barcodes:
            barcodes = all_barcodes
        else:
            barcodes = all_barcodeless

        for (project, library, lane, barcode, clust_count,
             clust_count_pf, nb_bases, nb_bases_r1_q30, nb_bases_r2_q30) in barcodes:
            reads_per_lane[lane] += clust_count_pf
            # For the moment, assume that nb_bases for r1 and r2 are the same.
            # TODO: remove this assumption by parsing ConversionStats.xml
            if not self.dataset.has_barcodes:
                barcode_info = self.barcodes_info.get('%s_%s' % (self.dataset.name, lane))
            else:
                barcode_info = self.barcodes_info.get('%s_%s_%s' % (self.dataset.name, lane, barcode))

            barcode_info[ELEMENT_NB_READS_SEQUENCED] = clust_count
            barcode_info[ELEMENT_NB_READS_PASS_FILTER] = clust_count_pf
            barcode_info[ELEMENT_NB_BASE_R1] = nb_bases
            barcode_info[ELEMENT_NB_BASE_R2] = nb_bases
            barcode_info[ELEMENT_NB_Q30_R1] = nb_bases_r1_q30
            barcode_info[ELEMENT_NB_Q30_R2] = nb_bases_r2_q30
        for run_element_id in self.barcodes_info:
            barcode = self.barcodes_info[run_element_id]
            reads_for_lane = reads_per_lane.get(barcode[ELEMENT_LANE])
            if reads_for_lane > 0:
                barcode[ELEMENT_PC_READ_IN_LANE] = barcode[ELEMENT_NB_READS_PASS_FILTER] / reads_for_lane

        for lane, barcode, clust_count in top_unknown_barcodes:
            unknown_element_id = '%s_%s_%s' % (self.dataset.name, lane, barcode)
            self.unexpected_barcodes[unknown_element_id] = {
                ELEMENT_RUN_ELEMENT_ID: unknown_element_id,
                ELEMENT_RUN_NAME: self.dataset.name,
                ELEMENT_LANE: lane,
                ELEMENT_PC_READ_IN_LANE: int(clust_count) / reads_per_lane.get(lane),
                ELEMENT_BARCODE: barcode,
                ELEMENT_NB_READS_PASS_FILTER: int(clust_count)
            }

    def send_data(self):
        return all(
            (
                pp('run_elements', self.barcodes_info.values(), ELEMENT_RUN_ELEMENT_ID),
                pp('unexpected_barcodes', self.unexpected_barcodes.values(), ELEMENT_RUN_ELEMENT_ID),
                pp('lanes', self.lanes.values(), ELEMENT_LANE_ID),
                pp('runs', [self.run], ELEMENT_RUN_NAME),
                pp('samples', self.libraries.values(), ELEMENT_SAMPLE_INTERNAL_ID, ['run_elements']),
                pp('projects', self.projects.values(), ELEMENT_PROJECT_ID, ['samples'])
            )
        )