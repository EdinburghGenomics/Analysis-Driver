from collections import Counter, defaultdict
import glob
import json
import os
from analysis_driver.clarity import get_sex_from_lims
from re import sub
from analysis_driver.app_logging import AppLogger
from analysis_driver.reader import demultiplexing_parsers, mapping_stats_parsers
from analysis_driver.reader.mapping_stats_parsers import parse_genotype_concordance
from analysis_driver.report_generation import rest_communication
from analysis_driver.config import default as cfg
from analysis_driver.report_generation import ELEMENT_RUN_NAME, ELEMENT_NUMBER_LANE, ELEMENT_RUN_ELEMENTS, \
    ELEMENT_BARCODE, ELEMENT_RUN_ELEMENT_ID, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_LIBRARY_INTERNAL_ID, \
    ELEMENT_LANE, ELEMENT_SAMPLES, ELEMENT_NB_READS_SEQUENCED, ELEMENT_NB_READS_PASS_FILTER, ELEMENT_NB_BASE_R1, \
    ELEMENT_NB_BASE_R2, ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2, ELEMENT_PC_READ_IN_LANE, ELEMENT_LANE_ID, \
    ELEMENT_PROJECT_ID, ELEMENT_SAMPLE_EXTERNAL_ID, ELEMENT_NB_READS_IN_BAM, ELEMENT_NB_MAPPED_READS, \
    ELEMENT_NB_DUPLICATE_READS, ELEMENT_NB_PROPERLY_MAPPED, ELEMENT_MEDIAN_COVERAGE, ELEMENT_PC_BASES_CALLABLE, \
    ELEMENT_LANE_NUMBER, ELEMENT_CALLED_GENDER, ELEMENT_PROVIDED_GENDER, ELEMENT_FREEMIX

gender_aliases = {'female': ['f', 'female', 'girl', 'women'],
                  'male': ['m', 'male', 'boy', 'man']}


def gender_alias(gender):
    g = str(gender).lower()
    for key in gender_aliases:
        if g in gender_aliases[key]:
            return key
    return 'unknown'


class Crawler(AppLogger):
    pass


class RunCrawler(Crawler):
    def __init__(self, run_id, samplesheet, conversion_xml_file=None):
        self.run_id = run_id
        self._populate_barcode_info_from_sample_sheet(samplesheet)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)

    def _populate_barcode_info_from_sample_sheet(self, samplesheet):
        self.barcodes_info = {}
        self.unexpected_barcodes = {}
        self.libraries = defaultdict(dict)
        self.lanes = defaultdict(dict)
        self.run = {ELEMENT_RUN_NAME: self.run_id, ELEMENT_RUN_ELEMENTS: []}
        self.projects = defaultdict(dict)

        for project_id, proj_obj in samplesheet.sample_projects.items():
            for sample_id_obj in proj_obj.sample_ids.values():
                for sample in sample_id_obj.samples:
                    for lane in sample.lane.split('+'):
                        barcode_info = {
                            ELEMENT_BARCODE: sample.barcode,
                            ELEMENT_RUN_ELEMENT_ID: '%s_%s_%s' % (self.run_id, lane, sample.barcode),
                            ELEMENT_RUN_NAME: self.run_id,
                            ELEMENT_PROJECT_ID: project_id,
                            ELEMENT_SAMPLE_INTERNAL_ID: sample.sample_id,
                            ELEMENT_LIBRARY_INTERNAL_ID: sample.sample_name,
                            ELEMENT_LANE: lane
                        }
                        self.barcodes_info[barcode_info[ELEMENT_RUN_ELEMENT_ID]] = barcode_info

                        # Populate the libraries
                        lib = self.libraries[sample.sample_name]
                        lib[ELEMENT_SAMPLE_INTERNAL_ID] = sample.sample_id
                        lib[ELEMENT_PROJECT_ID] = project_id
                        lib[ELEMENT_LIBRARY_INTERNAL_ID] = sample.sample_name
                        if ELEMENT_RUN_ELEMENTS not in lib:
                            lib[ELEMENT_RUN_ELEMENTS] = []
                        lib[ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        # Populate the projects
                        self.projects[project_id][ELEMENT_PROJECT_ID] = project_id
                        if ELEMENT_SAMPLES not in self.projects[project_id]:
                            self.projects[project_id][ELEMENT_SAMPLES] = set()
                        self.projects[project_id][ELEMENT_SAMPLES].add(sample.sample_id)

                        # Populate the lanes
                        lane_id = '%s_%s' % (self.run_id, lane)
                        ln = self.lanes[lane_id]
                        ln[ELEMENT_RUN_NAME] = self.run_id
                        ln[ELEMENT_LANE_ID] = lane_id
                        ln[ELEMENT_LANE_NUMBER] = int(lane)
                        if ELEMENT_RUN_ELEMENTS not in self.lanes[lane_id]:
                            ln[ELEMENT_RUN_ELEMENTS] = []
                        ln[ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        # Populate the run
                        self.run[ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        barcode_info = {
                            ELEMENT_BARCODE: 'unknown',
                            ELEMENT_RUN_ELEMENT_ID: '%s_%s_%s' % (self.run_id, lane, 'unknown'),
                            ELEMENT_RUN_NAME: self.run_id,
                            ELEMENT_PROJECT_ID: 'default',
                            ELEMENT_SAMPLE_INTERNAL_ID: 'Undetermined',
                            ELEMENT_LIBRARY_INTERNAL_ID: 'Undetermined',
                            ELEMENT_LANE: lane
                        }
                        self.barcodes_info[barcode_info[ELEMENT_RUN_ELEMENT_ID]] = barcode_info

        for project_id in self.projects:
            self.projects[project_id][ELEMENT_SAMPLES] = list(self.projects[project_id][ELEMENT_SAMPLES])

        # Add the unknown to the lane
        for lane_id in self.lanes:
            lane = self.lanes[lane_id][ELEMENT_LANE_NUMBER]
            unknown = '%s_%s_%s' % (self.run_id, lane, 'unknown')
            self.lanes[lane_id][ELEMENT_RUN_ELEMENTS].append(unknown)

        self.run[ELEMENT_NUMBER_LANE] = len(self.lanes)

    def _populate_barcode_info_from_conversion_file(self, conversion_xml_file):
        all_barcodes_per_lanes, top_unknown_barcodes_per_lanes = demultiplexing_parsers.parse_conversion_stats(conversion_xml_file)
        nb_read_per_lane = Counter()
        for project, library, lane, barcode, clust_count, clust_count_pf, nb_bases, nb_bases_r1_q30, nb_bases_r2_q30 in all_barcodes_per_lanes:
            nb_read_per_lane[lane] += int(clust_count_pf)
            barcode_info = self.barcodes_info.get('%s_%s_%s' % (self.run_id, lane, barcode))
            barcode_info[ELEMENT_NB_READS_SEQUENCED] = int(clust_count)
            barcode_info[ELEMENT_NB_READS_PASS_FILTER] = int(clust_count_pf)
            # For the paired end reads. For the moment, assume that r1 and r2 are the same lengths.
            # TODO: we should remove this assumption by parsing ConversionStats.xml
            barcode_info[ELEMENT_NB_BASE_R1] = int(nb_bases)
            barcode_info[ELEMENT_NB_BASE_R2] = int(nb_bases)
            barcode_info[ELEMENT_NB_Q30_R1] = int(nb_bases_r1_q30)
            barcode_info[ELEMENT_NB_Q30_R2] = int(nb_bases_r2_q30)
        for run_element_id in self.barcodes_info:
            run_element = self.barcodes_info[run_element_id]
            if nb_read_per_lane.get(run_element[ELEMENT_LANE]) > 0:
                run_element[ELEMENT_PC_READ_IN_LANE] = run_element[ELEMENT_NB_READS_PASS_FILTER] / nb_read_per_lane.get(run_element[ELEMENT_LANE])

        for lane, barcode, clust_count in top_unknown_barcodes_per_lanes:
            barcode_info = {
                ELEMENT_RUN_ELEMENT_ID: '%s_%s_%s' % (self.run_id, lane, barcode),
                ELEMENT_RUN_NAME: self.run_id,
                ELEMENT_LANE: lane,
                ELEMENT_PC_READ_IN_LANE: int(clust_count) / nb_read_per_lane.get(lane),
                ELEMENT_BARCODE: barcode,
                ELEMENT_NB_READS_PASS_FILTER: int(clust_count)
            }
            self.unexpected_barcodes[barcode_info[ELEMENT_RUN_ELEMENT_ID]] = barcode_info

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
        if not cfg.get('rest_api'):
            self.warn('rest_api is not set in the config: Cancel upload')
            return None
        
        return all(
            (
                rest_communication.post_or_patch('run_elements', self.barcodes_info.values(), ELEMENT_RUN_ELEMENT_ID),
                rest_communication.post_or_patch('unexpected_barcodes', self.unexpected_barcodes.values(), ELEMENT_RUN_ELEMENT_ID),
                rest_communication.post_or_patch('lanes', self.lanes.values(), ELEMENT_LANE_ID),
                rest_communication.post_or_patch('runs', [self.run], ELEMENT_RUN_NAME),
                rest_communication.post_or_patch('samples', self.libraries.values(), ELEMENT_SAMPLE_INTERNAL_ID, update_lists=['run_elements']),
                rest_communication.post_or_patch('projects', self.projects.values(), ELEMENT_PROJECT_ID, update_lists=['samples'])
            )
        )

class SampleCrawler(Crawler):
    def __init__(self, sample_id,  project_id, sample_dir):
        self.sample_id = sample_id
        self.project_id = project_id
        self.all_info = []
        self.sample = self._populate_lib_info(sample_dir)

    def _populate_lib_info(self, sample_dir):
        fastq_file = glob.glob(os.path.join(sample_dir, '*_R1.fastq.gz'))[0]
        external_sample_name = os.path.basename(fastq_file)[:-len('_R1.fastq.gz')]

        sample = {
            ELEMENT_SAMPLE_INTERNAL_ID: self.sample_id,
            ELEMENT_PROJECT_ID: self.project_id,
            ELEMENT_SAMPLE_EXTERNAL_ID: external_sample_name
        }

        bamtools_path = glob.glob(os.path.join(sample_dir, 'bamtools_stats.txt'))
        if bamtools_path:
            total_reads, mapped_reads, duplicate_reads, proper_pairs = mapping_stats_parsers.parse_bamtools_stats(bamtools_path[0])
            sample[ELEMENT_NB_READS_IN_BAM] = int(total_reads)
            sample[ELEMENT_NB_MAPPED_READS] = int(mapped_reads)
            sample[ELEMENT_NB_DUPLICATE_READS] = int(duplicate_reads)
            sample[ELEMENT_NB_PROPERLY_MAPPED] = int(proper_pairs)
        else:
            self.critical('Missing bamtools_stats.txt')

        yaml_metric_paths = glob.glob(
            os.path.join(
                sample_dir,
                '*%s-sort-highdepth-stats.yaml' % external_sample_name
            )
        )
        if yaml_metric_paths:
            yaml_metric_path = yaml_metric_paths[0]
            median_coverage = mapping_stats_parsers.parse_highdepth_yaml_file(yaml_metric_path)
            sample[ELEMENT_MEDIAN_COVERAGE] = median_coverage
        else:
            self.critical('Missing %s-sort-highdepth-stats.yaml' % external_sample_name)

        bed_file_paths = glob.glob(os.path.join(sample_dir, '*%s-sort-callable.bed' % external_sample_name))
        if bed_file_paths:
            bed_file_path = bed_file_paths[0]
            coverage_per_type = mapping_stats_parsers.parse_callable_bed_file(bed_file_path)
            callable_bases = coverage_per_type.get('CALLABLE')
            total = sum(coverage_per_type.values())
            sample[ELEMENT_PC_BASES_CALLABLE] = callable_bases/total
        else:
            self.critical('Missing *%s-sort-callable.bed' % external_sample_name)
        sex_file_paths = glob.glob(os.path.join(sample_dir, '%s.sex' % external_sample_name))
        if not sex_file_paths:
            sex_file_paths = glob.glob(os.path.join(sample_dir, '.qc', '%s.sex' % external_sample_name))
        if sex_file_paths:
            with open(sex_file_paths[0]) as open_file:
                gender = open_file.read().strip()
                gender_from_lims = get_sex_from_lims(self.sample_id)
                sample[ELEMENT_PROVIDED_GENDER] = gender_alias(gender_from_lims)
                sample[ELEMENT_CALLED_GENDER] = gender_alias(gender)

        self_sm_file = glob.glob(os.path.join(sample_dir, external_sample_name + '.selfSM'))
        if self_sm_file:
            with open(self_sm_file) as f:
                header = sub(r'  +', ' ', f.readline()).split(' ')
                data = sub('  +', ' ', f.readline()).split(' ')
                if len(header) != len(data):
                    self.warn('Inconsistent header and data length: %s and %s' % (len(header), len(data)))
                assert header[6] == 'FREEMIX'
                sample[ELEMENT_FREEMIX] = data[6]
        else:
            self.error('Missing selfSM file')

        genotype_validation_paths = glob.glob(os.path.join(sample_dir, '%s_genotype_validation.txt' % external_sample_name))
        if genotype_validation_paths:
            genotyping_results = parse_genotype_concordance(genotype_validation_paths[0])
            genotyping_result = genotyping_results.get(self.sample_id)
            if genotyping_result:
                sample.update(genotyping_result)
            else:
                self.critical('Sample %s not found in file %s' % (self.sample_id, genotype_validation_paths[0]))

        return sample

    def send_data(self):
        if not cfg.get('rest_api'):
            self.warn('rest_api is not set in the config: Cancel upload')
            return

        return rest_communication.post_or_patch('samples', [self.sample], ELEMENT_SAMPLE_INTERNAL_ID)
