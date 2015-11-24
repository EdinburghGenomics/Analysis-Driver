#!/usr/bin/env python
from collections import Counter, defaultdict
import json
import os
from analysis_driver.reader.demultiplexing_parsers import parse_conversion_stats
from analysis_driver.report_generation.rest_communication import post_entry, patch_entry
from analysis_driver.config import default as cfg

__author__ = 'tcezard'

ELEMENT_RUN_ELEMENT_ID = 'run_element_id'
ELEMENT_RUN_NAME = 'run_id'
ELEMENT_LANE = 'lane'
ELEMENT_BARCODE = 'barcode'
ELEMENT_PROJECT='project'
ELEMENT_LIBRARY_INTERNAL_ID = 'library_id'
ELEMENT_SAMPLE_INTERNAL_ID = 'sample_id'
ELEMENT_NB_READS_SEQUENCED = 'total_reads'
ELEMENT_NB_READS_PASS_FILTER = 'passing_filter_reads'
ELEMENT_PC_READ_IN_LANE = 'pc_reads_in_lane'
ELEMENT_NB_BASE_R1 = 'bases_r1'
ELEMENT_NB_BASE_R2 = 'bases_r2'
ELEMENT_NB_Q30_R1 = 'q30_bases_r1'
ELEMENT_NB_Q30_R2 = 'q30_bases_r2'
ELEMENT_FASTQC_REPORT_R1 = 'fastqc_report_r1'
ELEMENT_FASTQC_REPORT_R2 = 'fastqc_report_r2'

ELEMENT_SAMPLE_EXTERNAL_ID = 'user_sample_id'
ELEMENT_NB_READS_IN_BAM = 'bam_file_reads'
ELEMENT_NB_MAPPED_READS = 'mapped_reads'
ELEMENT_NB_SEC_MAPPED_READS = 'Nb secondary alignments'
ELEMENT_NB_DUPLICATE_READS = 'duplicate_reads'
ELEMENT_NB_PROPERLY_MAPPED = 'properly_mapped_reads'
ELEMENT_MEDIAN_INSERT_SIZE = 'Median Insert Size'
ELEMENT_MEAN_INSERT_SIZE = 'Mean Insert Size'
ELEMENT_MEDIAN_COVERAGE = 'median_coverage'
ELEMENT_PC_BASES_CALLABLE = 'pc_callable'
ELEMENT_MEAN_COVERAGE = 'Mean coverage'
ELEMENT_RUN_ELEMENTS = 'run_elements'

class RunCrawler:

    def __init__(self, run_id, samplesheet, conversion_xml_file=None):
        self.run_id = run_id
        self._populate_barcode_info_from_SampleSheet(samplesheet)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)

    def _populate_barcode_info_from_SampleSheet(self, samplesheet):
        self.barcodes_info={}
        self.libraries = defaultdict(dict)
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
                        if not ELEMENT_RUN_ELEMENTS in self.libraries[sample.sample_name]:
                            self.libraries[sample.sample_name][ELEMENT_RUN_ELEMENTS] = []
                        self.libraries[sample.sample_name][ELEMENT_RUN_ELEMENTS].append(barcode_info[ELEMENT_RUN_ELEMENT_ID])

                        barcode_info = {}
                        barcode_info[ELEMENT_BARCODE]='unknown'
                        barcode_info[ELEMENT_RUN_ELEMENT_ID] = '%s_%s_%s'%(self.run_id, lane, 'unknown')
                        barcode_info[ELEMENT_RUN_NAME]=self.run_id
                        barcode_info[ELEMENT_PROJECT] = 'default'
                        barcode_info[ELEMENT_SAMPLE_INTERNAL_ID]='Undetermined'
                        barcode_info[ELEMENT_LIBRARY_INTERNAL_ID]='Undetermined'
                        barcode_info[ELEMENT_LANE]=lane
                        self.barcodes_info[barcode_info[ELEMENT_RUN_ELEMENT_ID]]=(barcode_info)
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
            run_element[ELEMENT_PC_READ_IN_LANE] = run_element[ELEMENT_NB_READS_PASS_FILTER]  / \
                                                          nb_read_per_lane.get(run_element[ELEMENT_LANE])

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
            'demultiplexing' : self.barcodes_info.values(),
            'unexpected_barcodes' : self.unexpected_barcode_info.values(),
            'libraries' : self.libraries.values()
        }
        with open(json_file, 'w') as open_file:
            json.dump(payload, open_file)

    def write_json_per_sample(self, sample_dir):
        self.libraries.values()
        for library in self.libraries:
            file_name = os.path.join(sample_dir,self.libraries[library][ELEMENT_SAMPLE_INTERNAL_ID])
            if os.path.exists(file_name):
                with open(file_name) as open_file:
                    payload = json.load(open_file)
            else:
                payload = {}
            for run_element_id in self.libraries[library][ELEMENT_RUN_ELEMENTS]:
                payload[run_element_id] = self.barcodes_info[run_element_id][ELEMENT_NB_READS_SEQUENCED]
            with open(file_name, '') as open_file:
                json.dump(payload, open_file)

    def send_data(self):
        #Send run elements
        if cfg['rest_api']:
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
                    id = payload.pop(ELEMENT_RUN_ELEMENT_ID.key)
                    patch_entry(url, payload, **{ELEMENT_RUN_ELEMENT_ID.key:id})

            #Send samples information
            array_json = self.libraries.values()
            url=cfg.query('rest_api','url') + 'samples/'
            for payload in array_json:
                lib_id = {ELEMENT_LIBRARY_INTERNAL_ID:payload.get(ELEMENT_LIBRARY_INTERNAL_ID)}
                if not post_entry(url, payload):
                    patch_entry(url, payload, **lib_id)



    def __str__(self):
        return self.write_report_wiki()
