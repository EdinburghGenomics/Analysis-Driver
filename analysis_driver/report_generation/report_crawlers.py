from collections import Counter, defaultdict
import json
from analysis_driver import util
from analysis_driver.clarity import get_sex_from_lims, get_user_sample_name
from analysis_driver.app_logging import AppLogger
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.reader import demultiplexing_parsers, mapping_stats_parsers
from analysis_driver.reader.demultiplexing_parsers import get_fastqscreen_results
from analysis_driver.rest_communication import post_or_patch as pp
from analysis_driver.config import default as cfg
from analysis_driver.constants import ELEMENT_RUN_NAME, ELEMENT_NUMBER_LANE, ELEMENT_RUN_ELEMENTS, \
    ELEMENT_BARCODE, ELEMENT_RUN_ELEMENT_ID, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_LIBRARY_INTERNAL_ID, \
    ELEMENT_LANE, ELEMENT_SAMPLES, ELEMENT_NB_READS_SEQUENCED, ELEMENT_NB_READS_PASS_FILTER, ELEMENT_NB_BASE_R1, \
    ELEMENT_NB_BASE_R2, ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2, ELEMENT_PC_READ_IN_LANE, ELEMENT_LANE_ID, \
    ELEMENT_PROJECT_ID, ELEMENT_SAMPLE_EXTERNAL_ID, ELEMENT_NB_READS_IN_BAM, ELEMENT_NB_MAPPED_READS, \
    ELEMENT_NB_DUPLICATE_READS, ELEMENT_NB_PROPERLY_MAPPED, ELEMENT_MEDIAN_COVERAGE, ELEMENT_PC_BASES_CALLABLE, \
    ELEMENT_LANE_NUMBER, ELEMENT_CALLED_GENDER, ELEMENT_PROVIDED_GENDER, ELEMENT_NB_READS_CLEANED, ELEMENT_NB_Q30_R1_CLEANED, \
    ELEMENT_NB_BASE_R2_CLEANED, ELEMENT_NB_Q30_R2_CLEANED, ELEMENT_NB_BASE_R1_CLEANED, ELEMENT_GENOTYPE_VALIDATION, \
    ELEMENT_SPECIES_CONTAMINATION, ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL, ELEMENT_PCNT_UNMAPPED



class Crawler(AppLogger):
    def _check_config(self):
        if cfg.get('rest_api') and cfg.query('rest_api', 'url'):
            return True
        else:
            self.warn('rest_api is not configured. Cancel upload')
            return False


class RunCrawler(Crawler):
    def __init__(self, run_id, samplesheet, conversion_xml_file=None, run_dir=None):
        self.run_id = run_id
        self._populate_barcode_info_from_sample_sheet(samplesheet)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)
        if run_dir:
            self._populate_barcode_info_from_seqtk_fqchk_files(run_dir)

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

    def _populate_barcode_info_from_sample_sheet(self, samplesheet):
        self.barcodes_info = defaultdict(dict)
        self.unexpected_barcodes = {}
        self.libraries = defaultdict(dict)
        self.lanes = defaultdict(dict)
        self.run = {ELEMENT_RUN_NAME: self.run_id, ELEMENT_RUN_ELEMENTS: []}
        self.projects = defaultdict(dict)

        for project_id, proj_obj in samplesheet.sample_projects.items():
            for sample_id_obj in proj_obj.sample_ids.values():
                for sample in sample_id_obj.samples:
                    for lane in sample.lane.split('+'):
                        run_element_id = '%s_%s_%s' % (self.run_id, lane, sample.barcode)

                        self.barcodes_info[run_element_id] = {
                            ELEMENT_BARCODE: sample.barcode,
                            ELEMENT_RUN_ELEMENT_ID: run_element_id,
                            ELEMENT_RUN_NAME: self.run_id,
                            ELEMENT_PROJECT_ID: project_id,
                            ELEMENT_SAMPLE_INTERNAL_ID: sample.sample_id,
                            ELEMENT_LIBRARY_INTERNAL_ID: sample.sample_name,
                            ELEMENT_LANE: lane
                        }

                        # Populate the libraries
                        lib = self.libraries[sample.sample_name]
                        lib[ELEMENT_SAMPLE_INTERNAL_ID] = sample.sample_id
                        lib[ELEMENT_PROJECT_ID] = project_id
                        lib[ELEMENT_LIBRARY_INTERNAL_ID] = sample.sample_name
                        self._update_doc_list(lib, k=ELEMENT_RUN_ELEMENTS, v=run_element_id)

                        # Populate the projects
                        proj = self.projects[project_id]
                        proj[ELEMENT_PROJECT_ID] = project_id
                        self._update_doc_set(proj, k=ELEMENT_SAMPLES, v=sample.sample_id)

                        # Populate the lanes
                        lane_id = '%s_%s' % (self.run_id, lane)
                        ln = self.lanes[lane_id]
                        ln[ELEMENT_RUN_NAME] = self.run_id
                        ln[ELEMENT_LANE_ID] = lane_id
                        ln[ELEMENT_LANE_NUMBER] = int(lane)
                        self._update_doc_list(ln, k=ELEMENT_RUN_ELEMENTS, v=run_element_id)

                        # Populate the run
                        self.run[ELEMENT_RUN_ELEMENTS].append(run_element_id)

                        unknown_element_id = '%s_%s_%s' % (self.run_id, lane, 'unknown')
                        self.barcodes_info[unknown_element_id] = {
                            ELEMENT_BARCODE: 'unknown',
                            ELEMENT_RUN_ELEMENT_ID: unknown_element_id,
                            ELEMENT_RUN_NAME: self.run_id,
                            ELEMENT_PROJECT_ID: 'default',
                            ELEMENT_SAMPLE_INTERNAL_ID: 'Undetermined',
                            ELEMENT_LIBRARY_INTERNAL_ID: 'Undetermined',
                            ELEMENT_LANE: lane
                        }

        for project_id in self.projects:
            self.projects[project_id][ELEMENT_SAMPLES] = list(self.projects[project_id][ELEMENT_SAMPLES])

        # Add the unknown to the lane
        for lane_id in self.lanes:
            lane = self.lanes[lane_id][ELEMENT_LANE_NUMBER]
            unknown = '%s_%s_%s' % (self.run_id, lane, 'unknown')
            self.lanes[lane_id][ELEMENT_RUN_ELEMENTS].append(unknown)

        self.run[ELEMENT_NUMBER_LANE] = len(self.lanes)

    def _populate_barcode_info_from_seqtk_fqchk_files(self, run_dir):
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            if barcode_info[ELEMENT_BARCODE] == 'unknown':
                fq_chk_files = util.find_files(
                    run_dir,
                    'fastq',
                    'Undetermined_S0_L00%s_R*_001.fastq.gz.fqchk' % barcode_info[ELEMENT_LANE]
                )
            else:
                fq_chk_files = util.find_files(
                    run_dir,
                    'fastq',
                    barcode_info[ELEMENT_PROJECT_ID],
                    barcode_info[ELEMENT_SAMPLE_INTERNAL_ID],
                    '*_S*_L00%s_R*_001.fastq.gz.fqchk' % barcode_info[ELEMENT_LANE]
                )

            if len(fq_chk_files) == 2:
                fq_chk_files.sort()
                nb_read, nb_base, lo_q, hi_q = demultiplexing_parsers.parse_seqtk_fqchk_file(fq_chk_files[0], q_threshold=30)
                barcode_info[ELEMENT_NB_READS_CLEANED] = nb_read
                barcode_info[ELEMENT_NB_BASE_R1_CLEANED] = nb_base
                barcode_info[ELEMENT_NB_Q30_R1_CLEANED] = hi_q
                nb_read, nb_base, lo_q, hi_q = demultiplexing_parsers.parse_seqtk_fqchk_file(fq_chk_files[1], q_threshold=30)
                barcode_info[ELEMENT_NB_BASE_R2_CLEANED] = nb_base
                barcode_info[ELEMENT_NB_Q30_R2_CLEANED] = hi_q

            elif len(fq_chk_files) == 1:
                raise AnalysisDriverError('Only one fqchk file found in %s for %s' % (run_dir, run_element_id))

            elif barcode_info[ELEMENT_NB_READS_PASS_FILTER] == 0:
                barcode_info[ELEMENT_NB_READS_CLEANED] = 0
                barcode_info[ELEMENT_NB_BASE_R1_CLEANED] = 0
                barcode_info[ELEMENT_NB_Q30_R1_CLEANED] = 0
                barcode_info[ELEMENT_NB_BASE_R2_CLEANED] = 0
                barcode_info[ELEMENT_NB_Q30_R2_CLEANED] = 0

            else:
                raise AnalysisDriverError('%s fqchk files found in %s for %s' % (len(fq_chk_files), run_dir, run_element_id))

    def _populate_barcode_info_from_conversion_file(self, conversion_xml):
        all_barcodes, top_unknown_barcodes = demultiplexing_parsers.parse_conversion_stats(conversion_xml)
        reads_per_lane = Counter()

        for (project, library, lane, barcode, clust_count,
             clust_count_pf, nb_bases, nb_bases_r1_q30, nb_bases_r2_q30) in all_barcodes:
            reads_per_lane[lane] += clust_count_pf
            # For the moment, assume that nb_bases for r1 and r2 are the same.
            # TODO: remove this assumption by parsing ConversionStats.xml
            barcode_info = self.barcodes_info.get('%s_%s_%s' % (self.run_id, lane, barcode))
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
            unknown_element_id = '%s_%s_%s' % (self.run_id, lane, barcode)
            self.unexpected_barcodes[unknown_element_id] = {
                ELEMENT_RUN_ELEMENT_ID: unknown_element_id,
                ELEMENT_RUN_NAME: self.run_id,
                ELEMENT_LANE: lane,
                ELEMENT_PC_READ_IN_LANE: int(clust_count) / reads_per_lane.get(lane),
                ELEMENT_BARCODE: barcode,
                ELEMENT_NB_READS_PASS_FILTER: int(clust_count)
            }

    def write_json(self, json_file):
        payload = {
            'run_elements': list(self.barcodes_info.values()),
            'unexpected_barcodes': list(self.unexpected_barcodes.values()),
            'lanes': list(self.lanes.values()),
            'runs': self.run,
            'samples': list(self.libraries.values()),
            'project': list(self.projects.values())
        }
        with open(json_file, 'w') as open_file:
            json.dump(payload, open_file, indent=4)

    def send_data(self):
        if self._check_config():
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


class SampleCrawler(Crawler):
    gender_aliases = {'female': ['f', 'female', 'girl', 'woman'], 'male': ['m', 'male', 'boy', 'man']}

    def __init__(self, sample_id,  project_id,  sample_dir):
        self.sample_id = sample_id
        self.project_id = project_id
        self.all_info = []
        self.sample = self._populate_lib_info(sample_dir)

    @classmethod
    def _gender_alias(cls, gender):
        for key in cls.gender_aliases:
            if str(gender).lower() in cls.gender_aliases[key]:
                return key
        return 'unknown'

    def _populate_lib_info(self, sample_dir):
        external_sample_name = get_user_sample_name(self.sample_id, lenient=True)

        sample = {
            ELEMENT_SAMPLE_INTERNAL_ID: self.sample_id,
            ELEMENT_PROJECT_ID: self.project_id,
            ELEMENT_SAMPLE_EXTERNAL_ID: external_sample_name
        }

        bamtools_path = util.find_file(sample_dir, 'bamtools_stats.txt')
        if bamtools_path:
            (total_reads, mapped_reads,
             duplicate_reads, proper_pairs) = mapping_stats_parsers.parse_bamtools_stats(bamtools_path)

            sample[ELEMENT_NB_READS_IN_BAM] = total_reads
            sample[ELEMENT_NB_MAPPED_READS] = mapped_reads
            sample[ELEMENT_NB_DUPLICATE_READS] = duplicate_reads
            sample[ELEMENT_NB_PROPERLY_MAPPED] = proper_pairs
        else:
            self.critical('Missing bamtools_stats.txt')

        yaml_metric_path = util.find_file(sample_dir, '*%s-sort-highdepth-stats.yaml' % external_sample_name)
        if yaml_metric_path:
            median_coverage = mapping_stats_parsers.parse_highdepth_yaml_file(yaml_metric_path)
            sample[ELEMENT_MEDIAN_COVERAGE] = median_coverage
        else:
            self.critical('Missing %s-sort-highdepth-stats.yaml' % external_sample_name)

        bed_file_path = util.find_file(sample_dir, '*%s-sort-callable.bed' % external_sample_name)
        if bed_file_path:
            coverage_per_type = mapping_stats_parsers.parse_callable_bed_file(bed_file_path)
            callable_bases = coverage_per_type.get('CALLABLE')
            total = sum(coverage_per_type.values())
            sample[ELEMENT_PC_BASES_CALLABLE] = callable_bases/total
        else:
            self.critical('Missing *%s-sort-callable.bed' % external_sample_name)

        sex_file = util.find_file(sample_dir, '%s.sex' % external_sample_name)
        if not sex_file:
            sex_file = util.find_file(sample_dir, '.qc', '%s.sex' % external_sample_name)
        if sex_file:
            with open(sex_file) as f:
                gender = f.read().strip()
                gender_from_lims = get_sex_from_lims(self.sample_id)
                sample[ELEMENT_PROVIDED_GENDER] = self._gender_alias(gender_from_lims)
                sample[ELEMENT_CALLED_GENDER] = self._gender_alias(gender)

        genotype_validation_paths = util.find_files(sample_dir, '%s_genotype_validation.txt' % external_sample_name)
        if genotype_validation_paths:
            genotyping_results = mapping_stats_parsers.parse_genotype_concordance(genotype_validation_paths[0])
            genotyping_result = genotyping_results.get(self.sample_id)
            if genotyping_result:
                sample[ELEMENT_GENOTYPE_VALIDATION] = genotyping_result
            else:
                self.critical('Sample %s not found in file %s' % (self.sample_id, genotype_validation_paths[0]))


        species_contamination_paths = util.find_files(sample_dir, '%s_screen.txt' % external_sample_name)
        if species_contamination_paths:
            species_contamination_result = get_fastqscreen_results(species_contamination_paths[0], self.sample_id)
            if species_contamination_result:
                sample[ELEMENT_SPECIES_CONTAMINATION] = species_contamination_result
            else:
                self.critical('Contamination check unavailable for %s' % (self.sample_id))

        return sample

    def send_data(self):
        if self._check_config():
            return pp('samples', [self.sample], ELEMENT_SAMPLE_INTERNAL_ID)
