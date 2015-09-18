__author__ = 'mwham'
import csv
import os.path
from analysis_driver.app_logging import AppLogger, get_logger
from .run_info import RunInfo
from analysis_driver.config import default as cfg


app_logger = get_logger('reader')


def transform_sample_sheet(data_dir):
    before, header = _read_sample_sheet(data_dir, 'SampleSheet.csv')
    cols = before.readline().strip().split(',')
    after = open(os.path.join(data_dir, 'SampleSheet_analysis_driver.csv'), 'w')
    transformations = cfg['sample_sheet'].get('transformations', [])
    for idx, col in enumerate(cols):
        if col in transformations:
            cols[idx] = transformations[col]
    for line in header:
        after.write(line)
    after.write('[Data],\n')
    after.write(','.join(cols) + '\n')
    for line in before:
        after.write(line)
    before.close()


def _read_sample_sheet(data_dir, name):
    f = open(os.path.join(data_dir, name), 'r')
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

    return None, None


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
        error = 0
        if not self.run_info.barcode_len:
            self.run_info.critical('No barcode found in RunInfo.xml')
            error = 1
        if self.check_barcodes() != self.run_info.barcode_len:
            self.error(
                'Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' %
                (self.check_barcodes(), self.run_info.barcode_len)
            )
            error = 1
        self.debug('Done')
        return error

    def get_samples(self, sample_project, sample_id):
        return self.sample_projects[sample_project].sample_ids[sample_id].samples

    def _populate(self):
        f, header = _read_sample_sheet(self.data_dir, 'SampleSheet_analysis_driver.csv')
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        counter = 0
        for line in reader:
            if any(line):
                counter += 1
                sample_project = line[self._get_column(cols, 'sample_project')]
                sample_id = line[self._get_column(cols, 'sample_id')]

                new_sample = Sample(
                    sample_project=sample_project,
                    lane=line[self._get_column(cols, 'lane')],
                    sample_id=sample_id,
                    sample_name=line[self._get_column(cols, 'sample_name')],
                    barcode=line[self._get_column(cols, 'barcode')],
                    **self._get_all_cols(
                        line,
                        ignore=['sample_project', 'sample_id', 'sample_name' 'lane', 'barcode']
                    )
                )

                sample_project_obj = self._get_sample_project(sample_project)
                sample_id_obj = sample_project_obj.get_sample_id(sample_id)
                sample_id_obj.add_sample(new_sample)
        f.close()
        self.debug('Added %s samples' % counter)

    @staticmethod
    def _get_column(header, name):
        possible_fields = cfg['sample_sheet']['column_names'][name]
        for f in possible_fields:
            if f in header:
                return f
        return None

    @staticmethod
    def _get_all_cols(line_dict, ignore=None):
        d = {}
        for k, v in line_dict.items():
            if k not in ignore:
                d[k] = v
        return d

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
    """
    def __init__(self, sample_project, sample_id, sample_name, lane, barcode, **kwargs):
        self.sample_project = sample_project
        self.sample_id = sample_id
        self.sample_name = sample_name
        self.lane = lane
        self.barcode = barcode
        self.extra_data = kwargs
