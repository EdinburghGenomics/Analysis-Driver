import math
from itertools import islice
from collections import Counter, defaultdict
from xml.etree import ElementTree
from egcg_core import constants as c
from egcg_core.app_logging import logging_default as log_cfg

app_logger = log_cfg.get_logger(__name__)


def parse_json_stats(json_data):
    """Creates an array of tuples of run elements, unknown run elements and associated metadata."""
    all_run_elements = []
    unknown_run_elements = []

    for lane in json_data['ConversionResults']:
        for sample in lane['DemuxResults']:
            """ If the run is not multiplexed, the index_sequence will not be present in the JSON file, and the cluster
            count needs to be retrieved from a different variable in the JSON file. The index_sequence variable will be 
            prepopulated with 'barcodeless'. """
            index_sequence = "barcodeless"
            cluster_count_raw = lane['TotalClustersRaw']
            cluster_count_pf = lane['TotalClustersPF']
            if 'IndexMetrics' in sample:
                index_sequence = sample['IndexMetrics'][0]['IndexSequence']
                cluster_count_raw = sample['NumberReads']
                cluster_count_pf = sample['NumberReads']  # cluster count is equal to cluster count pf in multiplexed runs

            # obtaining number of bases r1 and r2 q30, confirming the read number first
            nb_bases_r1_q30, nb_bases_r2_q30 = 0, 0
            if sample['ReadMetrics'][0]['ReadNumber'] == 1 and sample['ReadMetrics'][1]['ReadNumber'] == 2:
                nb_bases_r1_q30 = sample['ReadMetrics'][0]['YieldQ30']
                nb_bases_r2_q30 = sample['ReadMetrics'][1]['YieldQ30']
            elif sample['ReadMetrics'][0]['ReadNumber'] == 2 and sample['ReadMetrics'][1]['ReadNumber'] == 1:
                nb_bases_r1_q30 = sample['ReadMetrics'][1]['YieldQ30']
                nb_bases_r2_q30 = sample['ReadMetrics'][0]['YieldQ30']

            if sample['ReadMetrics'][0]['Yield'] != sample['ReadMetrics'][1]['Yield']:
                # yield is expected to be equal for r1 and r2 in multiplexed and barcodeless runs
                raise NotImplementedError()

            all_run_elements.append((
                sample['SampleName'],
                lane['LaneNumber'],
                index_sequence,  # barcode sequence
                cluster_count_raw,  # cluster count
                cluster_count_pf,
                sample['ReadMetrics'][0]['Yield'],  # yield is equal for r1 and r2 in multiplexed and barcodeless runs
                nb_bases_r1_q30,
                nb_bases_r2_q30
            ))

    # parsing the unknown barcodes into their own array
    for lane in json_data['UnknownBarcodes']:
        for index, barcode in enumerate(lane['Barcodes']):
            unknown_run_elements.append((
                str(lane['Lane']),
                barcode,
                str(lane['Barcodes'][barcode])
            ))
            if index == 9:
                break

    return all_run_elements, unknown_run_elements


def parse_conversion_stats(xml_file, has_barcode):
    tree = ElementTree.parse(xml_file).getroot()
    all_barcodes_per_lanes = []
    all_barcodeless = []

    if not has_barcode:
        for project in tree.iter('Project'):
            if project.get('name') == 'all':
                continue

            for sample in project.findall('Sample'):
                if sample.get('name') == 'all':
                    continue

                for barcode in sample.findall('Barcode'):
                    if barcode.get('name') == 'all':
                        for lane in barcode.findall('Lane'):
                            clust_count = 0
                            clust_count_pf = 0
                            nb_bases = 0
                            nb_bases_r1_q30 = 0
                            nb_bases_r2_q30 = 0
                            for tile in lane.findall('Tile'):
                                clust_count += int(tile.find('Raw').find('ClusterCount').text)
                                clust_count_pf += int(tile.find('Pf').find('ClusterCount').text)
                                for read in tile.find('Pf').findall('Read'):
                                    if read.get('number') == '1':
                                        nb_bases += int(read.find('Yield').text)
                                        nb_bases_r1_q30 += int(read.find('YieldQ30').text)
                                    if read.get('number') == '2':
                                        nb_bases_r2_q30 += int(read.find('YieldQ30').text)
                            all_barcodeless.append(
                                (
                                    project.get('name'),
                                    sample.get('name'),
                                    lane.get('number'),
                                    barcode.get('name'),
                                    clust_count,
                                    clust_count_pf,
                                    nb_bases,
                                    nb_bases_r1_q30,
                                    nb_bases_r2_q30
                                )
                            )
    elif has_barcode:
        for project in tree.iter('Project'):
            if project.get('name') == 'all':
                continue

            for sample in project.findall('Sample'):
                if sample.get('name') == 'all':
                    continue

                for barcode in sample.findall('Barcode'):

                    if barcode.get('name') == 'all':
                        continue

                    for lane in barcode.findall('Lane'):
                        barcode.get('name')
                        clust_count = 0
                        clust_count_pf = 0
                        nb_bases = 0
                        nb_bases_r1_q30 = 0
                        nb_bases_r2_q30 = 0
                        for tile in lane.findall('Tile'):
                            clust_count += int(tile.find('Raw').find('ClusterCount').text)
                            clust_count_pf += int(tile.find('Pf').find('ClusterCount').text)
                            for read in tile.find('Pf').findall('Read'):
                                if read.get('number') == '1':
                                    nb_bases += int(read.find('Yield').text)
                                    nb_bases_r1_q30 += int(read.find('YieldQ30').text)
                                if read.get('number') == '2':
                                    nb_bases_r2_q30 += int(read.find('YieldQ30').text)

                        all_barcodes_per_lanes.append(
                            (
                                project.get('name'),
                                sample.get('name'),
                                lane.get('number'),
                                barcode.get('name'),
                                clust_count,
                                clust_count_pf,
                                nb_bases,
                                nb_bases_r1_q30,
                                nb_bases_r2_q30
                            )
                        )

    top_unknown_barcodes_per_lanes = []
    for lane in tree.find('Flowcell').findall('Lane'):
        for unknown_barcode in lane.iter('Barcode'):
            top_unknown_barcodes_per_lanes.append(
                (lane.get('number'), unknown_barcode.get('sequence'), unknown_barcode.get('count'))
            )
    return all_barcodes_per_lanes, top_unknown_barcodes_per_lanes, all_barcodeless


def parse_seqtk_fqchk_file(fqchk_file, q_threshold):
    with open(fqchk_file) as open_file:
        open_file.readline()
        header = open_file.readline().split()
        all_cycles = open_file.readline().split()
        first_cycle = open_file.readline().split()
        nb_read = int(first_cycle[1])
        nb_base = int(all_cycles[1])
        lo_q = 0
        hi_q = 0
        for i, h in enumerate(header[9:]):
            # header are %Q2
            if int(h[2:]) < q_threshold:
                lo_q += int(all_cycles[9+i])
            else:
                hi_q += int(all_cycles[9+i])
        return nb_read, nb_base, lo_q, hi_q


def parse_fastqscreen_file(fastqscreen_file, focal_species):
    """
    Parse fastqscreen's output file
    :return dict: the maximum number of reads mapped uniquely (singly or multiple times) to a contaminant species
    :return float: % reads unmapped to focal species
    :return float: % reads with no hits to any of the genomes provided
    :return int: number of reads mapped in total
    """
    if focal_species is None:
        app_logger.warning('No species name available')
        return None

    uniquely_mapped = {}
    # set to 100% as default in case no focal species is available
    focal_species_pc_unmapped = float(100)
    species = []
    with open(fastqscreen_file) as f:
        lines = f.readlines()
        total_reads_mapped = int(((lines[0]).split(': ')[2]).rstrip('\n'))
        hit_no_genomes = float((lines[-1]).split(': ')[1])
        results = (lines[2:-2])
        for r in results:
            species_name = r.split('\t')[0].replace('_', ' ')
            species.append(species_name)

    if focal_species not in species:
        app_logger.warning('The focal species is not included in the contaminant database')

    for r in results:
        species_name = r.split('\t')[0].replace('_', ' ')
        results = r.split('\t')[1:12]

        nb_uniquely_mapped = int(r.split('\t')[4]) + int(r.split('\t')[6])
        uniquely_mapped[species_name] = nb_uniquely_mapped
        if species_name == focal_species:
            focal_species_pc_unmapped = float(results[2])

    uniquely_mapped = {k: v for k, v in uniquely_mapped.items() if v != 0}

    fastqscreen_qc = {
        c.ELEMENT_CONTAMINANT_UNIQUE_MAP: uniquely_mapped,
        c.ELEMENT_TOTAL_READS_MAPPED: total_reads_mapped,
        c.ELEMENT_PCNT_UNMAPPED_FOCAL: focal_species_pc_unmapped,
        c.ELEMENT_PCNT_UNMAPPED: hit_no_genomes
    }
    app_logger.debug('Extracted fastqscreen QC: %s', fastqscreen_qc)
    return fastqscreen_qc


def read_histogram_file(histogram_file):
    histograms = defaultdict(Counter)
    with open(histogram_file) as f:
        for line in f:
            sp_line = line.split()
            if 'PATCH' in sp_line[0]:
                # Ignore the PATCH chromosome as they are full of Ns
                continue
            else:
                histograms[sp_line[0]][int(sp_line[1])] += int(sp_line[2])
    return histograms


def collapse_histograms(histograms):
    res = Counter()
    for histogram in histograms.values():
        res.update(histogram)
    return res


def get_percentile(histogram, percentile):
    """
    Find a percentile from a histogram.
    :param dict histogram:
    :param int percentile: a scalar between 0 and 100.
    :return: A value for the given percentile
    """
    n_percentile = None
    sorted_set_of_keys = sorted(histogram.keys())
    sum_of_values = sum(histogram.values())
    percentile_index = sum_of_values * (percentile / 100)
    index_sum = 0
    if len(sorted_set_of_keys) < 1:
        n_percentile = 0
    else:
        for i in range(len(sorted_set_of_keys)):
            k = sorted_set_of_keys[i]
            index_sum += histogram.get(k)
            if index_sum > percentile_index:
                n_percentile = k
                break
            elif index_sum == percentile_index:
                n_percentile = (k + sorted_set_of_keys[i + 1])/2
                break
    return n_percentile


def calculate_mean(histogram):
    sum_of_depths = sum(k * v for k, v in histogram.items())
    number_of_depths = sum(histogram.values())
    return sum_of_depths/number_of_depths


def calculate_median(histogram):
    return get_percentile(histogram, 50)


def calculate_sd(histogram):
    mean_depth = int(calculate_mean(histogram))
    number_of_depths = 0
    sum_of_squared_differences = 0
    for depth, count in histogram.items():
        number_of_depths += count
        sd = ((depth - mean_depth) ** 2) * count
        sum_of_squared_differences += sd
    return math.sqrt(sum_of_squared_differences/number_of_depths)


def calculate_bases_at_coverage(histogram):
    bases_5x = sum(histogram[i] for i in histogram.keys() if i > 5)
    bases_15x = sum(histogram[i] for i in histogram.keys() if i > 15)
    bases_30x = sum(histogram[i] for i in histogram.keys() if i > 30)
    return bases_5x, bases_15x, bases_30x


def calculate_size_genome(histogram):
    return sum(histogram.values())


def calculate_evenness(histogram):
    """
    Calculation of evenness based on http://www.nature.com/jhg/journal/v61/n7/full/jhg201621a.html.
    R code extracted from the paper:
    C=round(mean(D))
    D2=D[D<=C]
    E=1-(length(D2)-sum(D2)/C)/length(D)
    Where D is a vector of number representing the coverage at every bases.
    """
    rounded_mean = round(calculate_mean(histogram))
    if rounded_mean > 0:
        low_half_hist = {k: v for k, v in histogram.items() if k <= rounded_mean}
        a = sum(low_half_hist.values())
        b = sum(k * v for k, v in low_half_hist.items()) / rounded_mean
        evenness = 1 - ((a - b) / sum(histogram.values()))
    else:
        evenness = 0
    return evenness


def get_coverage_statistics(histogram_file):
    # Read the histogram file keeping each chrom separated
    histograms = read_histogram_file(histogram_file)
    # Collapse all chroms into one hist
    histogram = collapse_histograms(histograms)
    # Calculate statistics
    coverage_mean = calculate_mean(histogram)
    coverage_median = calculate_median(histogram)
    coverage_sd = calculate_sd(histogram)
    coverage_percentiles = {c.ELEMENT_PERCENTILE_5: get_percentile(histogram, 5),
                            c.ELEMENT_PERCENTILE_25: get_percentile(histogram, 25),
                            c.ELEMENT_PERCENTILE_50: get_percentile(histogram, 50),
                            c.ELEMENT_PERCENTILE_75: get_percentile(histogram, 75),
                            c.ELEMENT_PERCENTILE_95: get_percentile(histogram, 95)}

    bases_5x, bases_15x, bases_30x = calculate_bases_at_coverage(histogram)
    bases_at_coverage = {
        c.ELEMENT_BASES_AT_5X: bases_5x, c.ELEMENT_BASES_AT_15X: bases_15x, c.ELEMENT_BASES_AT_30X: bases_30x
    }

    genome_size = calculate_size_genome(histogram)
    evenness = calculate_evenness(histogram)

    return coverage_mean, coverage_median, coverage_sd, coverage_percentiles, bases_at_coverage, genome_size, evenness


def get_coverage_y_chrom(histogram_file, chr_name='chrY'):
    # Read the histogram file keeping each chrom separated
    histograms = read_histogram_file(histogram_file)
    # Calculate the mean coverage of Y chromosome
    if chr_name in histograms:
        return calculate_mean(histograms.get('chrY'))


def parse_welldup_file(welldup_file):
    dup_per_lane = {}
    in_summary = False
    lane = None

    with open(welldup_file) as f:
        for line in f:
            if line.startswith('LaneSummary:'):
                lane = int(line.split()[1])
                in_summary = True
            elif in_summary and line.startswith('Level: 3'):
                pc_dup = line.split()[12].lstrip('(').rstrip(')')
                dup_per_lane[lane] = round(float(pc_dup) * 100, 3)
                in_summary = False
            elif line.startswith('Lane: '):
                in_summary = False
    return dup_per_lane


def parse_adapter_trim_file(adapter_trim_file, run_id):
    adapters_trimmed_by_id = {}
    with open(adapter_trim_file) as open_file:
        open_file = open_file.read()
        adapter_trim_info = open_file.split('\n\n')[0].split('\n')
        for line in islice(adapter_trim_info, 1, None):
            lane, read, project, sample_id, sample_name, sample_number, trimmed_bases, pc_bases = line.split()
            run_element_info = (run_id, sample_id, lane)
            if not adapters_trimmed_by_id.get(run_element_info):
                adapters_trimmed_by_id[run_element_info] = {}
            adapters_trimmed_by_id[run_element_info]['read_%s_trimmed_bases' % read] = int(trimmed_bases)
    return adapters_trimmed_by_id


def parse_fastq_filterer_stats(filterer_stats):
    """Parse the key value pairs stats file produced by fastq filterer"""
    content = {}
    with open(filterer_stats) as open_file:
        for line in open_file:
                sp_line = line.strip().split()
                if sp_line and len(sp_line) > 1:
                    content[sp_line[0]] = sp_line[1]
    return content


def parse_interop_summary(summary_file):
    def grab_sections(f):
        sections = defaultdict(list)
        in_section = False
        section_name = None

        for line in f:
            if line.startswith('#'):
                continue
            sp_line = [s.strip() for s in line.strip().split(',')]
            if sp_line == ['']:
                in_section = False
            elif not in_section:
                section_name = sp_line[0]
                in_section = True
            elif sp_line[0].startswith('Read'):
                section_name = sp_line[0]
            else:
                sections[section_name].append(sp_line)
        return sections

    def parse_read_section(lines):
        values_per_lane = {}
        for sp_line in lines[1:]:
            if sp_line[0].isdigit() and sp_line[1] == '-':
                values_per_lane[sp_line[0]] = {
                    'pc_clust_pf': float(sp_line[4].split()[0]),  # avg Cluster PF
                    'pc_clust_pf_stdev': float(sp_line[4].split()[2]),  # std dev Cluster PF
                    'phasing': float(sp_line[5].split('/')[0].strip()),  # Phasing
                    'prephasing': float(sp_line[5].split('/')[1].strip()),  # Prephasing
                    'pc_q30': float(sp_line[10]),  # '%>=Q30'
                    'pc_aligned': float(sp_line[13].split()[0]),  # avg %aligned'
                    'pc_aligned_stdev': float(sp_line[13].split()[2]),  # std dev %aligned '
                    'pc_error': float(sp_line[14].split()[0]),  # 'avg %error'
                    'pc_error_stdev': float(sp_line[14].split()[2]),  # 'std dev %error'
                    'intensity_c1': float(sp_line[18].split()[0]),  # ' avg Intensity C1'
                    'intensity_c1_stdev': float(sp_line[18].split()[2]),  # 'std dev Intensity C1'
                }
        return values_per_lane

    with open(summary_file) as open_file:
        sections = grab_sections(open_file)
    metrics_per_lane = defaultdict(dict)
    if 'Read 1' not in sections or ('Read 2' not in sections and 'Read 3' not in sections):
        return dict(metrics_per_lane)
    for lane, metrics in parse_read_section(sections.pop('Read 1')).items():
        for metric in metrics:
            metrics_per_lane[lane][metric + '_r1'] = metrics[metric]

    if 'Read 2' in sections:
        read2_metrics_per_lane = parse_read_section(sections.pop('Read 2'))
    else:
        read2_metrics_per_lane = parse_read_section(sections.pop('Read 3'))
    for lane, metrics in read2_metrics_per_lane.items():
        for metric in metrics:
            metrics_per_lane[lane][metric + '_r2'] = metrics[metric]
    return dict(metrics_per_lane)
