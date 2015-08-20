__author__ = 'mwham'
import csv
import os.path
from analysis_driver.app_logging import AppLogger
from .run_info import RunInfo
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg


class SampleSheet(AppLogger):
    """
    Represents an instance of SampleSheet.csv. It ignores all lines until '[Data]' is found in the first
    column, then starts reading the CSV.
    """
    def __init__(self, data_dir):
        """
        :param str data_dir:  A file path to the input_data folder containing SampleSheet.csv
        """
        self.data_dir = data_dir

        self.run_info = RunInfo(self.data_dir)
        self.sample_projects = {}  # {name: samples} {str: Sample}
        self._populate()
        self.debug('Sample project entries: ' + str(self.sample_projects))

    def check_barcodes(self):
        """
        For each sample project, check that all the DNA barcodes are the same length
        :return: The DNA barcode length for the sample sheet
        :rtype: int
        """
        last_sample = None
        for name, sample_project in self.sample_projects.items():
            self.debug('Checking sample project ' + name)
            for name2, sample_id in sample_project.sample_ids.items():
                self.debug('Checking sample id ' + name2)
                for sample in sample_id.samples:
                    try:
                        if len(sample.barcode) == len(last_sample.barcode):
                            pass
                        else:
                            raise ValueError(
                                'Odd barcode length for %s: %s (%s) in sample project \'%s\' ' % (
                                    sample.id, sample.barcode, len(sample.barcode), name
                                )
                            )
                    except AttributeError:
                        pass
                    finally:
                        last_sample = sample
        self.debug('Barcode check done. Barcode len: %s' % len(last_sample.barcode))
        return len(last_sample.barcode)

    def generate_mask(self):
        self.debug('Generating mask...')
        mask = self.run_info.mask
        barcode_len = self.check_barcodes()

        out = ['y' + str(mask.num_cycles(mask.upstream_read) - 1) + 'n']

        for i in mask.index_lengths:
            diff = i - barcode_len
            out.append('i' + str(barcode_len) + 'n' * diff)

        out.append('y' + str(mask.num_cycles(mask.downstream_read) - 1) + 'n')
        self.debug(out)
        return ','.join(out)

    def validate(self):
        self.debug('Validating...')
        if not self.run_info.barcode_len:
            self.run_info.critical('No barcode found in RunInfo.xml')
        if self.check_barcodes() != self.run_info.barcode_len:
            self.error(
                'Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' %
                (self.check_barcodes(), self.run_info.barcode_len)
            )
        self.debug('Done')

    def _detect_format(self):
        f = self._read_sample_sheet()
        header = f.readline().strip().split(',')
        header.sort()
        for name, form in cfg['sample_sheet_formats'].items():
            if header == self._get_column_names(form):
                self.info('Detected sample sheet format: ' + name)
                return ColumnSet(**cfg['sample_sheet_formats'][name])

        raise AnalysisDriverError('No valid format for sample sheet found in config')

    @staticmethod
    def _get_column_names(form):
        names = []

        for name in form.values():
            if type(name) is str:
                names.append(name)
            elif type(name) is list:
                names.extend(name)
        names.sort()
        return names

    @staticmethod
    def _get_other_cols(line, form):
        d = {}
        for col, val in line.items():
            if col in form.other_cols:
                d[col] = val
        return d

    def _populate(self):
        form = self._detect_format()
        f = self._read_sample_sheet()
        reader = csv.DictReader(f)
        counter = 0
        for line in reader:
            if any(line):
                counter += 1
                sample_project = line[form.sample_project]
                sample_id = line[form.sample_id]

                new_sample = Sample(
                    sample_project=sample_project,
                    lane=line[form.lane],
                    sample_id=sample_id,
                    barcode=line[form.barcode],
                    **self._get_other_cols(line, form)
                )

                sample_project_obj = self._get_sample_project(sample_project)
                sample_id_obj = sample_project_obj.get_sample_id(sample_id)
                sample_id_obj.add_sample(new_sample)
        f.close()
        self.debug('Added %s samples' % counter)

    def _read_sample_sheet(self):
        f = open(os.path.join(self.data_dir, 'SampleSheet.csv'), 'r')
        counter = 1
        while not next(f).startswith('[Data]'):
            pass  # do nothing until [Data]
            counter += 1

        self.debug('Starting reading sample sheet from line ' + str(counter))
        return f

    def _get_sample_project(self, name):
        sample_project = ValueError('Could not add sample project ' + name)
        try:
            sample_project = self.sample_projects[name]
        except KeyError:
            sample_project = SampleProject(name)
            self.sample_projects[name] = sample_project
        finally:
            return sample_project


class SampleProject:
    """
    Represents a sample project, and contains a list of SampleID objects.
    """
    def __init__(self, name):
        """
        :param name: The sample project name
        """
        self.name = name
        self.sample_ids = {}

    def get_sample_id(self, name):
        sample_id = ValueError('Could not add sample id ' + name)
        try:
            sample_id = self.sample_ids[name]
        except KeyError:
            sample_id = SampleID(name, self.name)
            self.sample_ids[name] = sample_id
        finally:
            return sample_id


class SampleID:
    def __init__(self, name, sample_project):
        self.name = name
        self.sample_project = sample_project
        self.samples = []

    def add_sample(self, sample):
        assert sample.sample_project == self.sample_project and sample.sample_id == self.name,\
            'Adding invalid sample project to ' + self.name + ': ' + sample.sample_project
        self.samples.append(sample)


class Sample:
    """
    This represents a Sample, i.e. a line in SampleSheet.csv below the '[Data]' marker. Supports dict-style
    attribute fetching, e.g. sample_object['lane']
    For practical purposes, this has the same functionality as ColumnSet
    """
    def __init__(self, sample_project, sample_id, lane, barcode, **kwargs):
        self.sample_project = sample_project
        self.sample_id = sample_id
        self.land = lane
        self.barcode = barcode
        self.extra_data = kwargs


class ColumnSet:
    def __init__(self, **kwargs):
        self.data = {}
        for k, v in kwargs.items():
            self.data[k] = v

    def __getattr__(self, attr):
        return self.data[attr]

