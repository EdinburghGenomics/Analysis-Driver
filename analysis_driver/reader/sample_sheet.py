import csv
import os.path
from collections import defaultdict
from egcg_core.app_logging import AppLogger, logging_default as log_cfg
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver import config

app_logger = log_cfg.get_logger('reader')


def transform_sample_sheet(data_dir, remove_barcode=False):
    """
    Read SampleSheet.csv, translate column names and write to SampleSheet_analysis_driver.csv
    :param str data_dir: Full path to a data directory containing SampleSheet.csv
    :param bool remove_barcode: Whether to remove barcodes from the sample sheet
    """
    before, header = _read_sample_sheet(os.path.join(data_dir, 'SampleSheet.csv'))
    cols = before.readline().strip().split(',')
    out_lines = []
    transformations = config.sample_sheet_config.get('transformations', [])
    for idx, col in enumerate(cols):
        if col in transformations:
            cols[idx] = transformations[col]
    for line in header:
        out_lines.append(line.strip())
    out_lines.append('[Data],')
    out_lines.append(','.join(cols))
    index_col = index2_col = -1
    if remove_barcode:
        if 'Index' in cols:
            index_col = cols.index('Index')
        if 'Index2' in cols:
            index2_col = cols.index('Index2')
    for line in before:
        if remove_barcode:
            sp_line = line.strip().split(',')
            if index_col >= 0:
                sp_line[index_col] = ''
            if index2_col >= 0:
                sp_line[index2_col] = ''
            out_lines.append(','.join(sp_line))
        else:
            out_lines.append(line.strip())
    before.close()
    out_lines.append('')
    with open(os.path.join(data_dir, 'SampleSheet_analysis_driver.csv'), 'w') as after:
        after.write('\n'.join(out_lines))


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
        if line.startswith('[Data]'):
            app_logger.debug('Starting reading sample sheet from line ' + str(counter))
            return f, header
        else:
            counter += 1
            header.append(line)
    f.close()
    return None, None


class SampleSheet(AppLogger):
    def __init__(self, filename):
        self.projects = {}  # {name: samples} {str: Sample}
        self.filename = filename
        self.lanes = defaultdict(list)
        self._populate()
        self.debug('Sample project entries: ' + str(self.projects))

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

    def validate(self, reads):
        """
        Ensure that the SampleSheet is consistent with itself and RunInfo
        :param .run_info.Reads reads: Reads object to check against
        """
        self.debug('Validating...')
        if self.has_barcodes != reads.has_barcodes:
            raise AnalysisDriverError(
                'Barcodedness mismatch: %s (sample_sheet) and %s (run_info)' % (
                    self.has_barcodes, reads.has_barcodes
                )
            )
        if self.has_barcodes and self.barcode_len != reads.barcode_len:
            raise AnalysisDriverError(
                'Barcode mismatch: %s (sample_sheet) and %s (run_info)' % (
                    self.barcode_len, reads.barcode_len
                )
            )
        return True

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
        colnames = sample_sheet_line.keys()

        self.project_id = self._get_column_value(sample_sheet_line, 'sample_project', colnames)
        self.sample_id = self._get_column_value(sample_sheet_line, 'sample_id', colnames)
        self.sample_name = self._get_column_value(sample_sheet_line, 'sample_name', colnames)
        self.lanes = self._get_column_value(sample_sheet_line, 'lane', colnames).split('+')
        self.barcode = self._get_column_value(sample_sheet_line, 'barcode', colnames)

        self.extra_data = sample_sheet_line  # everything else has been removed by dict.pop

    @staticmethod
    def _get_column_name(colname, cols):
        possible_fields = config.sample_sheet_config['column_names'][colname]
        for field in possible_fields:
            if field in cols:
                return field
        raise KeyError('Could not find column for ' + colname)

    @classmethod
    def _get_column_value(cls, line, colname, cols):
        return line.pop(cls._get_column_name(colname, cols))
