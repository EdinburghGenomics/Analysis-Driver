#!/usr/bin/env python
from collections import  Counter
import xml.etree.ElementTree as ET

__author__ = 'tcezard'


def parse_demultiplexing_stats(xml_file):
    """parse the demultiplexing_stats.xml to extract number of read for each barcodes"""
    tree = ET.parse(xml_file).getroot()
    all_elements = []
    for project in tree.iter('Project'):
        if project.get('name') == 'default': continue
        for sample in project.findall('Sample'):
            for barcode in sample.findall('Barcode'):
                if project.get('name') != 'all' and barcode.get('name') == 'all': continue
                for lane in barcode.findall('Lane'):
                    all_elements.append((project.get('name'), sample.get('name'), barcode.get('name'),
                                         lane.get('number'), lane.find('BarcodeCount').text))
    return all_elements


def parse_conversion_stats(xml_file):
    tree = ET.parse(xml_file).getroot()
    all_barcodes_per_lanes = []

    for project in tree.iter('Project'):
        if project.get('name') == 'all': continue
        for sample in project.findall('Sample'):
            if sample.get('name') == 'all': continue
            for barcode in sample.findall('Barcode'):
                if barcode.get('name') == 'all': continue
                for lane in barcode.findall('Lane'):
                    barcode.get('name')
                    clust_count = 0
                    clust_count_pf = 0
                    nb_bases = 0
                    nb_bases_r1q30 = 0
                    nb_bases_r2q30 = 0
                    for tile in lane.findall('Tile'):
                        clust_count += int(tile.find('Raw').find('ClusterCount').text)
                        clust_count_pf += int(tile.find('Pf').find('ClusterCount').text)
                        for read in tile.find('Pf').findall('Read'):
                            if read.get('number') == "1":
                                nb_bases += int(read.find('Yield').text)
                                nb_bases_r1q30 += int(read.find('YieldQ30').text)
                            if read.get('number') == "2":
                                nb_bases_r2q30 += int(read.find('YieldQ30').text)
                    all_barcodes_per_lanes.append(
                        (project.get('name'), sample.get('name'), lane.get('number'), barcode.get('name'),
                         clust_count, clust_count_pf, nb_bases, nb_bases_r1q30, nb_bases_r2q30))
    top_unknown_barcodes_per_lanes = []
    for lane in tree.find('Flowcell').findall('Lane'):
        for unknown_barcode in lane.iter('Barcode'):
            top_unknown_barcodes_per_lanes.append((lane.get('number'), unknown_barcode.get('sequence'), unknown_barcode.get('count')))
    return all_barcodes_per_lanes, top_unknown_barcodes_per_lanes

