import csv
import os.path
from collections import defaultdict

import re
from datetime import date, datetime

from egcg_core.app_logging import AppLogger, logging_default as log_cfg
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import sample_sheet_config

app_logger = log_cfg.get_logger('reader')


def transform_sample_sheet(data_dir, seqlab2=True, remove_barcode=False):
    """
    Read SampleSheet.csv, translate column names and write to SampleSheet_analysis_driver.csv
    :param str data_dir: Full path to a data directory containing SampleSheet.csv
    :param bool seqlab2: Whether to use seqlab2 or seqlab1 column transformations
    :param bool remove_barcode: Whether to remove barcodes from the sample sheet
    """
    sample_sheet_in, header = _read_sample_sheet(os.path.join(data_dir, 'SampleSheet.csv'))
    cols_in = sample_sheet_in.readline().strip().split(',')
    col_transformations = sample_sheet_config['seqlab2'] if seqlab2 else sample_sheet_config['seqlab1']
    cols_out = [col_transformations.get(c, c) for c in cols_in]

    out_lines = [line.strip() for line in header]
    out_lines.append(','.join(cols_out))

    for line in sample_sheet_in:
        if remove_barcode:
            sp_line = line.strip().split(',')
            sp_line[cols_out.index('Index')] = ''  # use the index of the 'Index' column!
            sp_line[cols_out.index('Index2')] = ''
            out_lines.append(','.join(sp_line))
        else:
            out_lines.append(line.strip())

    sample_sheet_in.close()
    out_lines.append('')
    with open(os.path.join(data_dir, 'SampleSheet_analysis_driver.csv'), 'w') as sample_sheet_out:
        sample_sheet_out.write('\n'.join(out_lines))


def generate_samplesheet_from_lims(run_id, filename, index1=False):
    from egcg_core import clarity
    instance = SampleSheet(filename=filename)

    def find_pooling_step_for_artifact(art, max_iteration=10, expected_pooling_step_name=None):
        nb_iteration = 0
        while len(art.input_artifact_list()) == 1:
            art = art.input_artifact_list()[0]
            if nb_iteration == max_iteration:
                raise ValueError('Cannot find pooling step after %s iteraction' % max_iteration)
            nb_iteration += 1
        if expected_pooling_step_name and art.parent_process.type.name != expected_pooling_step_name:
            raise ValueError(
                'Mismatching Step name: %s != %s' % (expected_pooling_step_name, art.parent_process.type.name)
            )
        return art.input_artifact_list()

    run_process = clarity.get_run(run_id)
    flowcell = set(run_process.parent_processes()).pop().output_containers()[0]
    for lane in flowcell.placements:
        if len(flowcell.placements[lane].reagent_labels) > 1:
            artifacts = find_pooling_step_for_artifact(flowcell.placements[lane],
                                                       expected_pooling_step_name='Create PDP Pool')
        else:
            artifacts = [flowcell.placements[lane]]
        for artifact in artifacts:
            assert len(artifact.samples) == 1
            assert len(artifact.reagent_labels) == 1
            sample = artifact.samples[0]
            reagent_label = artifact.reagent_labels[0]
            match = re.match('(\w{4})-(\w{4}) \(([ATCG]{8})-([ATCG]{8})\)', reagent_label)
            line = {
                'Sample_Project': sample.project.name,
                'Sample_ID': sample.name,
                'Sample_Name': artifact.name,
                'Lane': lane.split(':')[0],
            }
            if index1:
                line['Index'] = match.group(3)
            else:
                line['Index'] = ''
            new_line = Line(line)
            for l in new_line.lanes:
                instance.lanes[l].append(new_line)
            project = instance._get_project(new_line.project_id)
            sample_id = project.get_sample_id(new_line.sample_id)
            sample_id.add_line(new_line)
    instance._validate()
    return instance

def _read_sample_sheet(sample_sheet):
    """
    Scan down a sample sheet until a [Data] line, then return the file object for further reading
    :return: The sample sheet file object, and all lines above [Data]
    """
    f = open(sample_sheet, 'r')
    app_logger.debug('Opened ' + f.name)
    counter = 1
    header = []
    for line in f:
        header.append(line)
        counter += 1
        if line.startswith('[Data]'):
            app_logger.debug('Starting reading sample sheet from line ' + str(counter))
            return f, header

    f.close()
    return None, None

def _generate_sample_sheet_header():
    return [
        '[Header]', 'Date, ' + datetime.now().strftime('%d/%m/%Y'), 'Workflow, Generate FASTQ Only', '',
       '[Reads]', '151', '151', '',  '[Settings]' 'Adapter, AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
        'AdapterRead2, AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', '', '[Data]'
    ]


class SampleSheet(AppLogger):

    def __init__(self, filename):
        self.projects = {}  # {name: samples} {str: Sample}
        self.filename = filename
        self.lanes = defaultdict(list)
        if os.path.exists(self.filename):
            self._populate()
            self.debug('Sample project entries: ' + str(self.projects))
            self._validate()

    def _validate(self):
        self.barcode_len = self._check_barcodes()
        self.has_barcodes = bool(self.barcode_len)
        if not self.has_barcodes:
            self._validate_one_sample_per_lane()

    def _check_barcodes(self):
        """
        For each sample project, check that all the DNA barcodes are the same length
        :return: The DNA barcode length for the sample sheet
        """
        last_line = None

        for p in self.projects.values():
            for s in p.sample_ids.values():
                for l in s.lines:
                    if last_line and len(last_line.barcode) != len(l.barcode):
                        raise AnalysisDriverError(
                            'Unexpected barcode length for %s: %s in project %s' % (
                                l.sample_id, l.barcode, p.name
                            )
                        )
                    last_line = l

        self.debug('Barcode check done. Barcode len: %s', len(last_line.barcode))
        return len(last_line.barcode)

    def _validate_one_sample_per_lane(self):
        for k, v in self.lanes.items():
            if len(v) > 1:
                raise AnalysisDriverError('Barcodeless sample sheet has %s samples in lane %s' % (len(v), k))


    def _populate(self):
        f, header = _read_sample_sheet(self.filename)
        reader = csv.DictReader(f)
        counter = 0

        for line in reader:
            if any(line):
                counter += 1
                new_line = Line(line)

                for l in new_line.lanes:
                    self.lanes[l].append(new_line)
                project = self._get_project(new_line.project_id)
                sample_id = project.get_sample_id(new_line.sample_id)
                sample_id.add_line(new_line)
        f.close()
        self.debug('Added %s samples', counter)

    def _get_project(self, name):
        if name not in self.projects:
            self.projects[name] = ProjectID(name)
        return self.projects[name]


class ProjectID:
    """Represents a project, grouping SampleID objects"""
    def __init__(self, name):
        self.name = name
        self.sample_ids = {}

    def get_sample_id(self, name):
        if name not in self.sample_ids:
            self.sample_ids[name] = SampleID(name, project_id=self.name)
        return self.sample_ids[name]


class SampleID:
    """
    Represents a sample in the sample sheet, which can have multiple barcodes, and therefore multiple lines in
    the sample sheet.
    """
    def __init__(self, name, project_id):
        self.name = name
        self.project_id = project_id
        self.lines = []

    def add_line(self, line):
        assert line.project_id == self.project_id and line.sample_id == self.name,\
            'Adding invalid sample project to ' + self.name + ': ' + line.project_id
        self.lines.append(line)


class Line:
    """Represents a line in SampleSheet.csv below the '[Data]' marker."""
    def __init__(self, sample_sheet_line):
        self.project_id = sample_sheet_line.pop('Sample_Project')
        self.sample_id = sample_sheet_line.pop('Sample_ID')
        self.sample_name = sample_sheet_line.pop('Sample_Name')
        self.lanes = sample_sheet_line.pop('Lane').split('+')
        self.barcode = sample_sheet_line.pop('Index')
        self.extra_data = sample_sheet_line

