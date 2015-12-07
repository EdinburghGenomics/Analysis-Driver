#!/usr/bin/env python
from collections import Counter, defaultdict
import glob
import json
import os
from analysis_driver.app_logging import AppLogger
from analysis_driver.reader.demultiplexing_parsers import parse_conversion_stats
from analysis_driver.reader.mapping_stats_parsers import parse_bamtools_stats, parse_highdepth_yaml_file, \
    parse_callable_bed_file
from analysis_driver.report_generation import ELEMENT_RUN_NAME, ELEMENT_NUMBER_LANE, ELEMENT_RUN_ELEMENTS, \
    ELEMENT_BARCODE, ELEMENT_RUN_ELEMENT_ID, ELEMENT_PROJECT, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_LIBRARY_INTERNAL_ID, \
    ELEMENT_LANE, ELEMENT_SAMPLES, ELEMENT_NB_READS_SEQUENCED, ELEMENT_NB_READS_PASS_FILTER, ELEMENT_NB_BASE_R1, \
    ELEMENT_NB_BASE_R2, ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2, ELEMENT_PC_READ_IN_LANE, ELEMENT_LANE_ID, \
    ELEMENT_PROJECT_ID, ELEMENT_SAMPLE_EXTERNAL_ID, ELEMENT_NB_READS_IN_BAM, ELEMENT_NB_MAPPED_READS, \
    ELEMENT_NB_DUPLICATE_READS, ELEMENT_NB_PROPERLY_MAPPED, ELEMENT_MEDIAN_COVERAGE, ELEMENT_PC_BASES_CALLABLE, \
    ELEMENT_LANE_NUMBER
from analysis_driver.report_generation.rest_communication import post_entry, patch_entry
from analysis_driver.config import default as cfg

__author__ = 'tcezard'



class RunCrawler(AppLogger):

    def __init__(self, run_id, samplesheet, conversion_xml_file=None):
        self.run_id = run_id
        self._populate_barcode_info_from_SampleSheet(samplesheet)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)

    def _populate_barcode_info_from_SampleSheet(self, samplesheet):
        self.barcodes_info={}
        self.libraries = defaultdict(dict)
        self.lanes = defaultdict(dict)
        self.run={ELEMENT_RUN_NAME : self.run_id,
                  ELEMENT_RUN_ELEMENTS : []}
        self.projects=defaultdict(dict)
        for project_id, proj_obj in samplesheet.sample_projects.items():
            for sample_id_obj in proj_obj.sample_ids.values():
                for sample in sample_id_obj.samples:
                    for lane in sample.lane.split('+'):
                        barcode_info = {}
                        barcode_info[ELEMENT_BARCODE]=sample.barcode
                        barcode_info[ELEMENT_RUN_ELEMENT_ID] = '%s_%s_%s'%(self.run_id, lane, sample.barcode)
                        barcode_info[ELEMENT_RUN_NAME]=self.run_id
                        barcode_info[ELEMENT_PROJECT]=project_id
                        barcode_info[ELEMENT_SAMPLE_INTERNAL_ID]=sample.sample_id
                        barcode_info[ELEMENT_LIBRARY_INTERNAL_ID]=sample.sample_name
                        barcode_info[ELEMENT_LANE]=lane
                        self.barcodes_info[barcode_info[ELEMENT_RUN_ELEMENT_ID]]=(barcode_info)

                        #Populate the libraries
                        self.libraries[sample.sample_name][ELEMENT_SAMPLE_INTERNAL_ID] = sample.sample_id
                        self.libraries[sample.sample_name][ELEMENT_PROJECT] = project_id
                        self.libraries[sample.sample_name][ELEMENT_LIBRARY_INTERNAL_ID] = sample.sample_name
                        if not ELEMENT_RUN_ELEMENTS in self.libraries[sample.sample_name]:
                            self.libraries[sample.sample_name][ELEMENT_RUN_ELEMENTS] = []
                        self.libraries[sample.sample_name][ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        #Populate the projects
                        self.projects[project_id][ELEMENT_PROJECT_ID]=project_id
                        if not ELEMENT_SAMPLES in self.projects[project_id]:
                            self.projects[project_id][ELEMENT_SAMPLES] = []
                        self.projects[project_id][ELEMENT_SAMPLES].append(sample.sample_id)

                        #Populate the lanes
                        lane_id = '%s_%s'%(self.run_id, lane)
                        self.lanes[lane_id][ELEMENT_RUN_NAME] = self.run_id
                        self.lanes[lane_id][ELEMENT_LANE_ID] = lane_id
                        self.lanes[lane_id][ELEMENT_LANE_NUMBER] = int(lane)
                        if not ELEMENT_RUN_ELEMENTS in self.lanes[lane_id]:
                            self.lanes[lane_id][ELEMENT_RUN_ELEMENTS] = []
                        self.lanes[lane_id][ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        #Populate the run
                        self.run[ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        barcode_info = {}
                        barcode_info[ELEMENT_BARCODE]='unknown'
                        barcode_info[ELEMENT_RUN_ELEMENT_ID] = '%s_%s_%s'%(self.run_id, lane, 'unknown')
                        barcode_info[ELEMENT_RUN_NAME]=self.run_id
                        barcode_info[ELEMENT_PROJECT] = 'default'
                        barcode_info[ELEMENT_SAMPLE_INTERNAL_ID]='Undetermined'
                        barcode_info[ELEMENT_LIBRARY_INTERNAL_ID]='Undetermined'
                        barcode_info[ELEMENT_LANE]=lane
                        self.barcodes_info[barcode_info[ELEMENT_RUN_ELEMENT_ID]]=(barcode_info)
        self.run[ELEMENT_NUMBER_LANE] = len(self.lanes)
        self.unexpected_barcode_info={}

    def _populate_barcode_info_from_conversion_file(self, conversion_xml_file):
        all_barcodes_per_lanes, top_unknown_barcodes_per_lanes = parse_conversion_stats(conversion_xml_file)
        nb_read_per_lane=Counter()
        for project, library, lane, barcode, clust_count, clust_count_pf, nb_bases,\
            nb_bases_r1q30, nb_bases_r2q30, in all_barcodes_per_lanes:
            barcode_info = self.barcodes_info.get('%s_%s_%s'%(self.run_id, lane, barcode))
            barcode_info[ELEMENT_NB_READS_SEQUENCED]=int(clust_count)
            barcode_info[ELEMENT_NB_READS_PASS_FILTER]=int(clust_count_pf)
            nb_read_per_lane[lane]+=int(clust_count_pf)
            #For the paired end reads
            barcode_info[ELEMENT_NB_BASE_R1]=int(nb_bases)
            barcode_info[ELEMENT_NB_BASE_R2]=int(nb_bases)
            barcode_info[ELEMENT_NB_Q30_R1]=int(nb_bases_r1q30)
            barcode_info[ELEMENT_NB_Q30_R2]=int(nb_bases_r2q30)
        for run_element_id in self.barcodes_info:
            run_element = self.barcodes_info[run_element_id]
            if nb_read_per_lane.get(run_element[ELEMENT_LANE]) > 0 :
                run_element[ELEMENT_PC_READ_IN_LANE] = run_element[ELEMENT_NB_READS_PASS_FILTER] / nb_read_per_lane.get(run_element[ELEMENT_LANE])
            else:
                run_element[ELEMENT_PC_READ_IN_LANE] = ''

        for lane, barcode, clust_count in top_unknown_barcodes_per_lanes:
            barcode_info = {}
            barcode_info[ELEMENT_RUN_ELEMENT_ID] = '%s_%s_%s'%(self.run_id, lane, barcode)
            barcode_info[ELEMENT_RUN_NAME]=self.run_id
            barcode_info[ELEMENT_LANE]=lane
            barcode_info[ELEMENT_BARCODE]=barcode
            barcode_info[ELEMENT_NB_READS_PASS_FILTER]=int(clust_count)
            self.unexpected_barcode_info[barcode_info[ELEMENT_RUN_ELEMENT_ID]]=(barcode_info)

    def write_json(self, json_file):
        payload ={
            'run_elements' : list(self.barcodes_info.values()),
            'unexpected_barcodes' : list(self.unexpected_barcode_info.values()),
            'lanes' : list(self.lanes.values()),
            'runs' : self.run,
            'samples' : list(self.libraries.values()),
            'project' : list(self.projects.values())
        }
        with open(json_file, 'w') as open_file:
            json.dump(payload, open_file, indent=4)

    def update_json_per_sample(self, sample_dir):
        self.libraries.values()
        for library in self.libraries:
            file_name = os.path.join(sample_dir,self.libraries[library][ELEMENT_SAMPLE_INTERNAL_ID])
            if os.path.exists(file_name):
                with open(file_name) as open_file:
                    payload = json.load(open_file)
            else:
                payload = {}
            for run_element_id in self.libraries[library][ELEMENT_RUN_ELEMENTS]:
                payload[run_element_id] = self.barcodes_info[run_element_id]
            with open(file_name, 'w') as open_file:
                json.dump(payload, open_file, indent=4)

    def send_data(self):

        if not cfg.get('rest_api'):
            self.warn('rest_api is not set in the config: Cancel upload')
            return

        #Send run elements
        array_json = self.barcodes_info.values()
        url=cfg.query('rest_api','url') + 'run_elements/'
        for payload in array_json:
            if not post_entry(url, payload):
                id = payload.pop(ELEMENT_RUN_ELEMENT_ID)
                patch_entry(url, payload, **{ELEMENT_RUN_ELEMENT_ID:id})

        #Send unexpected barcodes
        array_json = self.unexpected_barcode_info.values()
        url=cfg.query('rest_api','url') + 'unexpected_barcodes/'
        for payload in array_json:
            if not post_entry(url, payload):
                id = payload.pop(ELEMENT_RUN_ELEMENT_ID)
                patch_entry(url, payload, **{ELEMENT_RUN_ELEMENT_ID:id})

        #Send lanes
        array_json = self.lanes.values()
        url=cfg.query('rest_api','url') + 'lanes/'
        for payload in array_json:
            if not post_entry(url, payload):
                id = payload.pop(ELEMENT_LANE_ID)
                patch_entry(url, payload, **{ELEMENT_LANE_ID:id})

        #Send runs
        url = cfg.query('rest_api', 'url') + 'runs/'
        payload = self.run
        if not post_entry(url, payload):
            patch_entry(url, payload)


        #Send samples information
        array_json = self.libraries.values()
        url=cfg.query('rest_api','url') + 'samples/'
        for payload in array_json:
            lib_id = {ELEMENT_LIBRARY_INTERNAL_ID:payload.get(ELEMENT_LIBRARY_INTERNAL_ID)}
            if not post_entry(url, payload):
                patch_entry(url, payload, **lib_id)

        #Send projects
        array_json = self.projects.values()
        url=cfg.query('rest_api','url') + 'projects/'
        for payload in array_json:
            if not post_entry(url, payload):
                id = payload.pop(ELEMENT_PROJECT_ID)
                patch_entry(url, payload, **{ELEMENT_PROJECT_ID:id})




class SampleCrawler(AppLogger):

    def __init__(self, sample_id,  project_id,  sample_dir):
        self.sample_id = sample_id
        self.project_id = project_id
        self.all_info = []
        self.sample = self._populate_lib_info(sample_dir)

    def _populate_lib_info(self, sample_dir):
        sample = {}
        sample[ELEMENT_SAMPLE_INTERNAL_ID]= self.sample_id
        sample[ELEMENT_PROJECT]= self.project_id
        fastq_file = glob.glob(os.path.join(sample_dir,"*_R1.fastq.gz"))[0]
        external_sample_name = os.path.basename(fastq_file)[:-len("_R1.fastq.gz")]
        sample[ELEMENT_SAMPLE_EXTERNAL_ID]= external_sample_name
        
        bamtools_path = glob.glob(os.path.join(sample_dir, 'bamtools_stats.txt'))
        if bamtools_path:
            total_reads, mapped_reads, duplicate_reads, proper_pairs = parse_bamtools_stats(bamtools_path[0])
            sample[ELEMENT_NB_READS_IN_BAM]= int(total_reads)
            sample[ELEMENT_NB_MAPPED_READS]= int(mapped_reads)
            sample[ELEMENT_NB_DUPLICATE_READS]= int(duplicate_reads)
            sample[ELEMENT_NB_PROPERLY_MAPPED]= int(proper_pairs)
        else:
            self.critical('Missing bamtools_stats.txt')

        yaml_metric_paths = glob.glob(os.path.join(sample_dir, '*%s-sort-highdepth-stats.yaml'%external_sample_name))
        if yaml_metric_paths:
            yaml_metric_path = yaml_metric_paths[0]
            median_coverage  = parse_highdepth_yaml_file(yaml_metric_path)
            sample[ELEMENT_MEDIAN_COVERAGE]= median_coverage
        else:
            self.critical('Missing %s-sort-highdepth-stats.yaml'%external_sample_name)

        bed_file_paths = glob.glob(os.path.join(sample_dir,'*%s-sort-callable.bed'%external_sample_name))
        if bed_file_paths:
            bed_file_path = bed_file_paths[0]
            coverage_per_type = parse_callable_bed_file(bed_file_path)
            callable_bases = coverage_per_type.get('CALLABLE')
            total = sum(coverage_per_type.values())
            sample[ELEMENT_PC_BASES_CALLABLE]= callable_bases/total
        else:
            self.critical('Missing *%s-sort-callable.bed'%external_sample_name)
        return sample

        self.run_id = run_id
        self._populate_barcode_info_from_SampleSheet(samplesheet)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)

    def write_json(self, json_file):
        payload ={
            'samples' : list(self.sample)
        }
        with open(json_file, 'w') as open_file:
            json.dump(payload, open_file, indent=4)


    def send_data(self):

        if not cfg.get('rest_api'):
            self.warn('rest_api is not set in the config: Cancel upload')
            return

        #Send sample
        payload = self.sample
        url=cfg.query('rest_api','url') + 'samples/'
        if not post_entry(url, payload):
            id = payload.get(ELEMENT_SAMPLE_INTERNAL_ID)
            patch_entry(url, payload, **{ELEMENT_SAMPLE_INTERNAL_ID:id})
