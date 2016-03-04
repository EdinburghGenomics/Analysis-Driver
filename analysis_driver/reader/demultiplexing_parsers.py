__author__ = 'tcezard'

import os
import sys
from xml.etree import ElementTree

sys.path.append('../..')
from analysis_driver.clarity import get_species_from_sample
from analysis_driver.constants import ELEMENT_SPECIES_CONTAMINATION, ELEMENT_CONTAMINANT_UNIQUE_MAP, ELEMENT_PCNT_UNMAPPED_FOCAL, ELEMENT_PCNT_UNMAPPED, ELEMENT_TOTAL_READS_MAPPED
from analysis_driver.app_logging import get_logger
app_logger = get_logger(__name__)



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


def parse_conversion_stats(xml_file):
    tree = ElementTree.parse(xml_file).getroot()
    all_barcodes_per_lanes = []

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
    return all_barcodes_per_lanes, top_unknown_barcodes_per_lanes


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


def get_focal_species(sample_id):
    myFocalSpecies = get_species_from_sample(sample_id)
    return myFocalSpecies

def parse_fastqscreen_file(filename, myFocalSpecies):
    """
    parse the fastq screen outfile
    :return int: the maximum number of reads mapped uniquely (singly or multiple times) to a contaminant species
    :return float: % reads unmapped to focal Species
    :return float: % reads with no hits to any of the genomes provided
    :return int: number of reads mapped in total
    """
    file = open(filename)
    lines = file.readlines()

    contaminantsUniquelyMapped = {}
    focalSpeciesPercentUnmapped = ''
    total_reads_mapped = int(((lines[0]).split(': ')[2]).rstrip('\n'))
    Hit_no_genomes = float((lines[-1]).split(': ')[1])
    speciesResults = (lines[2:-2])
    speciesList = []
    for result in speciesResults:
        speciesName = result.split('\t')[0]
        speciesName = speciesName.replace('_',' ')
        speciesList.append(speciesName)

    if myFocalSpecies in speciesList:
        for result in speciesResults:
            speciesName = result.split('\t')[0]
            speciesName = speciesName.replace('_',' ')
            speciesResults = result.split('\t')[1:12]

            if speciesName != myFocalSpecies:
                numberUniquelyMapped = int(result.split('\t')[4]) + int(result.split('\t')[6])
                contaminantsUniquelyMapped[speciesName] = numberUniquelyMapped
            elif speciesName == myFocalSpecies:
                focalSpeciesPercentUnmapped = float(speciesResults[2])
        contaminantsUniquelyMapped = {k:v for k,v in contaminantsUniquelyMapped.items() if v != 0}
        ELEMENT_SPECIES_CONTAMINATION = {ELEMENT_CONTAMINANT_UNIQUE_MAP:contaminantsUniquelyMapped, ELEMENT_TOTAL_READS_MAPPED:total_reads_mapped, ELEMENT_PCNT_UNMAPPED_FOCAL:focalSpeciesPercentUnmapped, ELEMENT_PCNT_UNMAPPED:Hit_no_genomes}
        return ELEMENT_SPECIES_CONTAMINATION
    else:
        app_logger.warning('The focal species is not included in the contaminant database')
        return [100, 100, 100]

def get_fastqscreen_results(filename, sample_id):
    myFocalSpecies = get_focal_species(sample_id)
    if myFocalSpecies == None:
        app_logger.warning('No species name available')
        return [100, 100, 100]
    else:
        fastqscreen_results = parse_fastqscreen_file(filename, myFocalSpecies)
        return fastqscreen_results

