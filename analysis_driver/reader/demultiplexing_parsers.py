import math
from collections import Counter, defaultdict
from xml.etree import ElementTree
from egcg_core.clarity import get_species_from_sample
from egcg_core.constants import ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL,\
    ELEMENT_PCNT_UNMAPPED, ELEMENT_TOTAL_READS_MAPPED

from egcg_core.app_logging import logging_default as log_cfg
app_logger = log_cfg.get_logger(__name__)


def parse_demultiplexing_stats(xml_file):
    """Parse the demultiplexing_stats.xml to extract number of read for each barcodes"""
    tree = ElementTree.parse(xml_file).getroot()
    all_elements = []
    for project in tree.iter('Project'):
        if project.get('name') == 'default':
            continue

        for sample in project.findall('Sample'):
            for barcode in sample.findall('Barcode'):
                if project.get('name') != 'all' and barcode.get('name') == 'all':
                    continue

                for lane in barcode.findall('Lane'):
                    all_elements.append(
                        (
                            project.get('name'),
                            sample.get('name'),
                            barcode.get('name'),
                            lane.get('number'),
                            lane.find('BarcodeCount').text
                        )
                    )
    return all_elements


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
                            barcode.get('name')
                            clust_count = 0
                            clust_count_pf = 0
                            nb_bases = 0
                            nb_bases_r1_q30 = 0
                            nb_bases_r2_q30 = 0
                            for tile in lane.findall('Tile'):
                                clust_count += int(tile.find('Raw').find('ClusterCount').text)
                                clust_count_pf += int(tile.find('Pf').find('ClusterCount').text)
                                # FIXME: Read numbers in the ConversionStats.xml are wrong when using barcodeless run
                                # Need to be fixed in bcl2fast and then changed here
                                for read in tile.find('Pf').findall('Read'):
                                    if read.get('number') == "2":
                                        nb_bases += int(read.find('Yield').text)
                                        nb_bases_r1_q30 += int(read.find('YieldQ30').text)
                                    if read.get('number') != "2":
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
                                if read.get('number') == "1":
                                    nb_bases += int(read.find('Yield').text)
                                    nb_bases_r1_q30 += int(read.find('YieldQ30').text)
                                if read.get('number') == "2":
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
        first_line = open_file.readline()
        header = open_file.readline().split()
        all_cycles = open_file.readline().split()
        first_cycle = open_file.readline().split()
        nb_read = int(first_cycle[1])
        nb_base = int(all_cycles[1])
        lo_q = 0
        hi_q = 0
        for i, h in enumerate(header[9:]):
            #header are %Q2
            if int(h[2:]) < q_threshold:
                lo_q += int(all_cycles[9+i])
            else:
                hi_q += int(all_cycles[9+i])
        return  nb_read, nb_base, lo_q, hi_q


def parse_fastqscreen_file(filename, myFocalSpecies):
    """
    parse the fastq screen outfile
    :return dict: the maximum number of reads mapped uniquely (singly or multiple times) to a contaminant species
    :return float: % reads unmapped to focal Species
    :return float: % reads with no hits to any of the genomes provided
    :return int: number of reads mapped in total
    """
    uniquelyMapped = {}
    focalSpeciesPercentUnmapped = ''
    speciesList = []
    with open(filename) as open_file:
        lines = open_file.readlines()
        total_reads_mapped = int(((lines[0]).split(': ')[2]).rstrip('\n'))
        Hit_no_genomes = float((lines[-1]).split(': ')[1])
        speciesResults = (lines[2:-2])
        for result in speciesResults:
            speciesName = result.split('\t')[0]
            speciesName = speciesName.replace('_',' ')
            speciesList.append(speciesName)

    if myFocalSpecies in speciesList:
        for result in speciesResults:
            speciesName = result.split('\t')[0]
            speciesName = speciesName.replace('_',' ')
            speciesResults = result.split('\t')[1:12]


            numberUniquelyMapped = int(result.split('\t')[4]) + int(result.split('\t')[6])
            uniquelyMapped[speciesName] = numberUniquelyMapped
            if speciesName == myFocalSpecies:
                focalSpeciesPercentUnmapped = float(speciesResults[2])
        uniquelyMapped = {k:v for k,v in uniquelyMapped.items() if v != 0}
        fastqscreen_result = {ELEMENT_CONTAMINANT_UNIQUE_MAP:uniquelyMapped,
                                         ELEMENT_TOTAL_READS_MAPPED:total_reads_mapped,
                                         ELEMENT_PCNT_UNMAPPED_FOCAL:focalSpeciesPercentUnmapped,
                                         ELEMENT_PCNT_UNMAPPED:Hit_no_genomes}
        return fastqscreen_result
    else:
        app_logger.warning('The focal species is not included in the contaminant database')
        fastqscreen_result = None
        return fastqscreen_result

def get_fastqscreen_results(filename, sample_id):
    myFocalSpecies = get_species_from_sample(sample_id)
    if myFocalSpecies is None:
        app_logger.warning('No species name available')
        return None
    else:
        fastqscreen_results = parse_fastqscreen_file(filename, myFocalSpecies)
        return fastqscreen_results

def read_histogram_file(input_file):
    histograms = defaultdict(Counter)
    with open(input_file) as openfile:
        for line in openfile:
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



def get_percentiles(histogram, percentile):
    """
    Calculate the percentiles of the histogram.
    :param percentiles: a scalar or an array of scalar comprised in between 0 and 100.
    :return A value or a array value for the given percentiles.
    """
    n_percentile=None
    sorted_set_of_keys = sorted(histogram.keys())
    sum_of_values = sum(histogram.values())
    percentile_index = sum_of_values * (percentile / 100)
    index_sum = 0
    if len(sorted_set_of_keys)<1:
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
    sumOfDepths = sum([k * v for k, v in histogram.items()])
    numberOfDepths = sum(histogram.values())
    return sumOfDepths/numberOfDepths

def calculate_median(histogram):
    return get_percentiles(histogram, 50)

def calculate_sd(histogram):
    meanDepth = int(calculate_mean(histogram))
    numberOfDepths = 0
    sumOfSquaredDifference = 0
    for depth, count  in histogram.items():
        numberOfDepths += count
        sd = ((depth - meanDepth) ** 2) * count
        sumOfSquaredDifference += sd
    return math.sqrt(sumOfSquaredDifference/numberOfDepths)

def get_coverage_statistics(histogram_file):
    # Read the histogram file keeping each chrom separated
    histograms = read_histogram_file(histogram_file)
    # Collapse all chroms into one hist
    histogram = collapse_histograms(histograms)
    # Calculate statistics
    coverage_mean = calculate_mean(histogram)
    coverage_median = calculate_median(histogram)
    coverage_sd = calculate_sd(histogram)
    return coverage_mean, coverage_median, coverage_sd

def get_coverage_Y_chrom(histogram_file, chr_name='chrY'):
    # Read the histogram file keeping each chrom separated
    histograms = read_histogram_file(histogram_file)
    # Calculate the mean coverage of Y chromosome
    if chr_name in histograms:
        return calculate_mean(histograms.get('chrY'))

def parse_welldup_file(welldup_file):
    dup_per_lane = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0}
    in_summary = False
    with open(welldup_file) as open_file:
        for line in open_file:
            if line.startswith('LaneSummary:'):
                lane = int(line.split()[1])
                in_summary=True
            elif in_summary == 1 and line.startswith('Level: 3'):
                pc_dup = line.split()[12].strip('(').rstrip(')')
                dup_per_lane[lane]=round(float(pc_dup)*100,3)
                in_summary = False
            elif line.startswith('Lane: '):
                in_summary = False
    return dup_per_lane



