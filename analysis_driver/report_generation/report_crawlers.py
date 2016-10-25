import json
from collections import Counter, defaultdict
from egcg_core import util
from egcg_core.app_logging import AppLogger
from egcg_core.rest_communication import post_or_patch as pp
from egcg_core import clarity
from analysis_driver.exceptions import PipelineError
from analysis_driver.reader import demultiplexing_parsers, mapping_stats_parsers
from analysis_driver.reader.demultiplexing_parsers import get_coverage_statistics, \
    parse_welldup_file, get_coverage_Y_chrom, parse_fastqscreen_file
from analysis_driver.reader.mapping_stats_parsers import parse_and_aggregate_genotype_concordance,\
    parse_vbi_selfSM, parse_vcf_stats
from analysis_driver.config import default as cfg
from egcg_core.constants import ELEMENT_RUN_NAME, ELEMENT_NUMBER_LANE, ELEMENT_RUN_ELEMENTS, \
    ELEMENT_BARCODE, ELEMENT_RUN_ELEMENT_ID, ELEMENT_SAMPLE_INTERNAL_ID, ELEMENT_LIBRARY_INTERNAL_ID, \
    ELEMENT_LANE, ELEMENT_SAMPLES, ELEMENT_NB_READS_SEQUENCED, ELEMENT_NB_READS_PASS_FILTER, ELEMENT_NB_BASE_R1, \
    ELEMENT_NB_BASE_R2, ELEMENT_NB_Q30_R1, ELEMENT_NB_Q30_R2, ELEMENT_PC_READ_IN_LANE, ELEMENT_LANE_ID, \
    ELEMENT_PROJECT_ID, ELEMENT_SAMPLE_EXTERNAL_ID, ELEMENT_NB_READS_IN_BAM, ELEMENT_NB_MAPPED_READS, \
    ELEMENT_NB_DUPLICATE_READS, ELEMENT_NB_PROPERLY_MAPPED, ELEMENT_MEDIAN_COVERAGE, ELEMENT_PC_BASES_CALLABLE, \
    ELEMENT_LANE_NUMBER, ELEMENT_CALLED_GENDER, ELEMENT_PROVIDED_GENDER, ELEMENT_NB_READS_CLEANED, ELEMENT_NB_Q30_R1_CLEANED, \
    ELEMENT_SPECIES_CONTAMINATION, ELEMENT_NB_BASE_R2_CLEANED, ELEMENT_NB_Q30_R2_CLEANED, ELEMENT_NB_BASE_R1_CLEANED, \
    ELEMENT_GENOTYPE_VALIDATION, ELEMENT_COVERAGE_STATISTICS, ELEMENT_MEAN_COVERAGE, ELEMENT_COVERAGE_PERCENTILES, \
    ELEMENT_BASES_AT_COVERAGE, ELEMENT_MEDIAN_COVERAGE_SAMTOOLS, ELEMENT_COVERAGE_SD, ELEMENT_FREEMIX, ELEMENT_SAMPLE_CONTAMINATION, \
    ELEMENT_GENDER_VALIDATION, ELEMENT_GENDER_HETX, ELEMENT_LANE_PC_OPT_DUP, ELEMENT_GENDER_COVY, ELEMENT_SNPS_TI_TV, \
    ELEMENT_SNPS_HET_HOM, ELEMENT_ADAPTER_TRIM_R1, ELEMENT_ADAPTER_TRIM_R2, ELEMENT_SAMPLE_PLATE, ELEMENT_SAMPLE_SPECIES, \
    ELEMENT_SAMPLE_EXPECTED_YIELD, ELEMENT_SAMPLE_EXPECTED_COVERAGE

_gender_aliases = {'female': ['f', 'female', 'girl', 'woman'], 'male': ['m', 'male', 'boy', 'man']}


def gender_alias(gender):
    for key in _gender_aliases:
        if str(gender).lower() in _gender_aliases[key]:
            return key
    return 'unknown'


def get_sample_information_from_lims(sample_name):
    lims_sample = clarity.get_sample(sample_name)
    gender = gender_alias(clarity.get_sample_gender(sample_name))
    plate_id, well = clarity.get_plate_id_and_well(sample_name)
    species = clarity.get_species_from_sample(sample_name)
    external_sample_name = clarity.get_user_sample_name(sample_name, lenient=True)
    yield_q30 = clarity.get_expected_yield_for_sample(sample_name)
    coverage = lims_sample.udf.get('Coverage', '')
    return {
        ELEMENT_SAMPLE_EXTERNAL_ID: external_sample_name,
        ELEMENT_SAMPLE_PLATE: plate_id,
        ELEMENT_PROVIDED_GENDER: gender,
        ELEMENT_SAMPLE_SPECIES: species,
        ELEMENT_SAMPLE_EXPECTED_YIELD: yield_q30,
        ELEMENT_SAMPLE_EXPECTED_COVERAGE: coverage
    }


class Crawler(AppLogger):
    def _check_config(self):
        if cfg.get('rest_api') and cfg.query('rest_api', 'url'):
            return True
        else:
            self.warning('rest_api is not configured. Cancel upload')
            return False

class RunCrawler(Crawler):
    def __init__(self, run_id, samplesheet, adapter_trim_file=None, conversion_xml_file=None, run_dir=None):
        self.run_id = run_id
        self.adapter_trim_file = adapter_trim_file
        self.samplesheet = samplesheet
        self._populate_barcode_info_from_sample_sheet(samplesheet)
        self._populate_from_lims()
        if adapter_trim_file:
            self._populate_barcode_info_from_adapter_file(adapter_trim_file)
        if conversion_xml_file:
            self._populate_barcode_info_from_conversion_file(conversion_xml_file)
        if run_dir:
            self._populate_barcode_info_from_seqtk_fqchk_files(run_dir)
        if run_dir:
            welldup_files = util.find_files(run_dir, run_id + '.wellduplicate')
            if welldup_files:
                self._populate_barcode_info_from_well_dup(welldup_files[0])

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
                        if not samplesheet.has_barcode:
                            run_element_id = '%s_%s' % (self.run_id, lane)
                            self.barcodes_info[run_element_id] = {
                                ELEMENT_RUN_ELEMENT_ID: run_element_id,
                                ELEMENT_RUN_NAME: self.run_id,
                                ELEMENT_PROJECT_ID: project_id,
                                ELEMENT_SAMPLE_INTERNAL_ID: sample.sample_id,
                                ELEMENT_LIBRARY_INTERNAL_ID: sample.sample_name,
                                ELEMENT_LANE: lane
                            }
                        else:
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

                        if samplesheet.has_barcode:
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

        if samplesheet.has_barcode:
            # Add the unknown to the lane
            for lane_id in self.lanes:
                lane = self.lanes[lane_id][ELEMENT_LANE_NUMBER]
                unknown = '%s_%s_%s' % (self.run_id, lane, 'unknown')
                self.lanes[lane_id][ELEMENT_RUN_ELEMENTS].append(unknown)

        self.run[ELEMENT_NUMBER_LANE] = len(self.lanes)

    def _populate_from_lims(self):
        for libname in self.libraries:
            self.libraries[libname].update(
                get_sample_information_from_lims(self.libraries[libname][ELEMENT_SAMPLE_INTERNAL_ID])
            )

    def _run_sample_lane_to_barcode(self, adapters_trimmed_by_id, has_barcode):
        run_element_adapters_trimmed = {}
        for adapter_id in adapters_trimmed_by_id:
            run_element_id = None
            run_id, sample_id, lane = adapter_id
            if has_barcode:
                for i in self.barcodes_info:
                    if self.barcodes_info[i][ELEMENT_RUN_NAME] == run_id \
                            and self.barcodes_info[i][ELEMENT_SAMPLE_INTERNAL_ID] == sample_id \
                            and self.barcodes_info[i][ELEMENT_LANE] == lane:
                        run_element_id = self.barcodes_info[i][ELEMENT_RUN_ELEMENT_ID]
                        break
                    else:
                        run_element_id = '%s_%s' % (run_id, lane)
            run_element_adapters_trimmed[run_element_id] = adapters_trimmed_by_id[adapter_id]

            return run_element_adapters_trimmed

    def _populate_barcode_info_from_adapter_file(self, adapter_trim_file):
        has_barcode = self.samplesheet.has_barcode
        parsed_trimmed_adapters = demultiplexing_parsers.parse_adapter_trim_file(adapter_trim_file, self.run_id)
        run_element_adapters_trimmed = self._run_sample_lane_to_barcode(parsed_trimmed_adapters, has_barcode)
        for run_element_id in run_element_adapters_trimmed:
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
                nb_read, nb_base, lo_q, hi_q = demultiplexing_parsers.parse_seqtk_fqchk_file(fq_chk_files[0], q_threshold=30)
                barcode_info[ELEMENT_NB_READS_CLEANED] = nb_read
                barcode_info[ELEMENT_NB_BASE_R1_CLEANED] = nb_base
                barcode_info[ELEMENT_NB_Q30_R1_CLEANED] = hi_q
                nb_read, nb_base, lo_q, hi_q = demultiplexing_parsers.parse_seqtk_fqchk_file(fq_chk_files[1], q_threshold=30)
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

    def _populate_barcode_info_from_well_dup(self, welldup_file):
        dup_per_lane = parse_welldup_file(welldup_file)
        for run_element_id in self.barcodes_info:
            barcode_info = self.barcodes_info.get(run_element_id)
            lane = barcode_info.get(ELEMENT_LANE)
            if int(lane) in dup_per_lane:
                barcode_info[ELEMENT_LANE_PC_OPT_DUP] = dup_per_lane.get(int(lane))

    def _populate_barcode_info_from_conversion_file(self, conversion_xml):
        all_barcodes, top_unknown_barcodes, all_barcodeless = demultiplexing_parsers.parse_conversion_stats(conversion_xml, self.samplesheet.has_barcode)
        reads_per_lane = Counter()
        if not self.samplesheet.has_barcode:
            barcodes = all_barcodeless
        else:
            barcodes = all_barcodes

        for (project, library, lane, barcode, clust_count,
             clust_count_pf, nb_bases, nb_bases_r1_q30, nb_bases_r2_q30) in barcodes:
            reads_per_lane[lane] += clust_count_pf
            # For the moment, assume that nb_bases for r1 and r2 are the same.
            # TODO: remove this assumption by parsing ConversionStats.xml
            if not self.samplesheet.has_barcode:
                barcode_info = self.barcodes_info.get('%s_%s' % (self.run_id, lane))
            else:
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

    def __init__(self, sample_id,  project_id,  sample_dir):
        self.sample_id = sample_id
        self.project_id = project_id
        self.all_info = []
        self.sample = self._populate_lib_info(sample_dir)

    @staticmethod
    def search_file(sample_dir, file_name):
        path_to_search = [
            [sample_dir],
            [sample_dir, '.qc']
        ]
        for path in path_to_search:
            path.append(file_name)
            f = util.find_file(*path)
            if f:
                return f

    def _populate_lib_info(self, sample_dir):

        sample = {
            ELEMENT_SAMPLE_INTERNAL_ID: self.sample_id,
            ELEMENT_PROJECT_ID: self.project_id,
        }

        sample.update(get_sample_information_from_lims(self.sample_id))
        external_sample_name = sample.get(ELEMENT_SAMPLE_EXTERNAL_ID)

        bamtools_path = self.search_file(sample_dir, 'bamtools_stats.txt')
        if bamtools_path:
            (total_reads, mapped_reads,
             duplicate_reads, proper_pairs) = mapping_stats_parsers.parse_bamtools_stats(bamtools_path)

            sample[ELEMENT_NB_READS_IN_BAM] = total_reads
            sample[ELEMENT_NB_MAPPED_READS] = mapped_reads
            sample[ELEMENT_NB_DUPLICATE_READS] = duplicate_reads
            sample[ELEMENT_NB_PROPERLY_MAPPED] = proper_pairs
        else:
            self.critical('Missing bamtools_stats.txt')

        samtools_path = self.search_file(sample_dir, 'samtools_stats.txt')
        if samtools_path:
            (total_reads, mapped_reads,
             duplicate_reads, proper_pairs) = mapping_stats_parsers.parse_samtools_stats(samtools_path)

            sample[ELEMENT_NB_READS_IN_BAM] = total_reads
            sample[ELEMENT_NB_MAPPED_READS] = mapped_reads
            sample[ELEMENT_NB_DUPLICATE_READS] = duplicate_reads
            sample[ELEMENT_NB_PROPERLY_MAPPED] = proper_pairs
        else:
            self.critical('Missing samtools_stats.txt')

        bed_file_path = self.search_file(sample_dir, '*%s-sort-callable.bed' % external_sample_name)
        if bed_file_path:
            coverage_per_type = mapping_stats_parsers.parse_callable_bed_file(bed_file_path)
            callable_bases = coverage_per_type.get('CALLABLE')
            total = sum(coverage_per_type.values())
            sample[ELEMENT_PC_BASES_CALLABLE] = callable_bases/total
        else:
            self.critical('Missing *%s-sort-callable.bed', external_sample_name)

        sex_file_path = self.search_file(sample_dir, '%s.sex' % external_sample_name)
        if sex_file_path:
            with open(sex_file_path) as f:
                gender, het_x = f.read().strip().split()
                sample[ELEMENT_CALLED_GENDER] = gender_alias(gender)
                sample[ELEMENT_GENDER_VALIDATION] = {ELEMENT_GENDER_HETX: het_x}

        genotype_validation_path = self.search_file(sample_dir, '%s_genotype_validation.txt' % external_sample_name)
        if genotype_validation_path:
            genotyping_results = parse_and_aggregate_genotype_concordance(genotype_validation_path)
            genotyping_result = genotyping_results.get(self.sample_id)
            if genotyping_result:
                sample[ELEMENT_GENOTYPE_VALIDATION] = genotyping_result
            else:
                self.critical('Sample %s not found in file %s', self.sample_id, genotype_validation_path)

        species_contamination_path = self.search_file(sample_dir, '%s_R1_screen.txt' % external_sample_name)
        if species_contamination_path:
            species_contamination_result = parse_fastqscreen_file(species_contamination_path,
                                                                   sample[ELEMENT_SAMPLE_SPECIES])
            if species_contamination_result:
                sample[ELEMENT_SPECIES_CONTAMINATION] = species_contamination_result
            else:
                self.critical('Contamination check unavailable for %s', self.sample_id)

        sample_contamination_path = self.search_file(sample_dir, '%s-chr22-vbi.selfSM' % external_sample_name)
        if not sample_contamination_path:
            sample_contamination_path = self.search_file(sample_dir, '%s-chr22-vbi.selfSM' % self.sample_id)
        if sample_contamination_path:
            freemix = parse_vbi_selfSM(sample_contamination_path)
            if freemix is not None:
                sample[ELEMENT_SAMPLE_CONTAMINATION] = {ELEMENT_FREEMIX: freemix}
            else:
                self.critical('freemix results from validateBamId are not available for %s', self.sample_id)

        coverage_statistics_path = self.search_file(sample_dir, '%s.depth' % external_sample_name)
        if coverage_statistics_path:
            mean, median, sd, coverage_percentiles, bases_at_coverage = get_coverage_statistics(coverage_statistics_path)
            coverage_statistics = {ELEMENT_MEAN_COVERAGE: mean,
                                   ELEMENT_MEDIAN_COVERAGE_SAMTOOLS: median,
                                   ELEMENT_COVERAGE_SD: sd,
                                   ELEMENT_COVERAGE_PERCENTILES: coverage_percentiles,
                                   ELEMENT_BASES_AT_COVERAGE: bases_at_coverage}
            sample[ELEMENT_COVERAGE_STATISTICS] = coverage_statistics
            sample[ELEMENT_MEDIAN_COVERAGE] = median
            if ELEMENT_GENDER_VALIDATION in sample:
                covY = get_coverage_Y_chrom(coverage_statistics_path)
                if covY: sample[ELEMENT_GENDER_VALIDATION][ELEMENT_GENDER_COVY] = covY
        else:
            self.critical('coverage statistics unavailable for %s', self.sample_id)

        vcf_stats_path = self.search_file(sample_dir, '%s.vcf.stats' % external_sample_name)
        if vcf_stats_path:
            ti_tv, het_hom = parse_vcf_stats(vcf_stats_path)
            if ELEMENT_SAMPLE_CONTAMINATION in sample:
                sample[ELEMENT_SAMPLE_CONTAMINATION][ELEMENT_SNPS_TI_TV] = ti_tv
                sample[ELEMENT_SAMPLE_CONTAMINATION][ELEMENT_SNPS_HET_HOM] = het_hom
            else:
                sample[ELEMENT_SAMPLE_CONTAMINATION]={ELEMENT_SNPS_TI_TV: ti_tv, ELEMENT_SNPS_HET_HOM:het_hom}

        return sample

    def send_data(self):
        if self._check_config():
            return pp('samples', [self.sample], ELEMENT_SAMPLE_INTERNAL_ID)