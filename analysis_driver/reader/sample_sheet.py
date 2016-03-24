from collections import defaultdict

from analysis_driver.exceptions import AnalysisDriverError

__author__ = 'mwham'
import csv
import os.path
from analysis_driver.app_logging import AppLogger, get_logger
from analysis_driver import config

app_logger = get_logger('reader')


def transform_sample_sheet(data_dir, remove_barcode=False):
    """
    Read SampleSheet.csv, translate column names and write to SampleSheet_analysis_driver.csv
    :param data_dir: Full path to a data directory containing SampleSheet.csv
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
    :return tuple[file, list] f, header: The sample sheet file object, and all lines above [Data]
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
    def __init__(self, filename, has_barcode=True):
        self.sample_projects = {}  # {name: samples} {str: Sample}
        self.filename = filename
        self._populate()
        self.debug('Sample project entries: ' + str(self.sample_projects))
        if not has_barcode:
            self.has_barcode = False
            self._validate_one_sample_per_lane()
        else:
            self.has_barcode = True


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


    def _validate_one_sample_per_lane(self):
        """
        Check that only one sample is present in each lane and raise AnalysisDriverError if more than one is found.
        """
        lane2samples = defaultdict(list)
        for name, sample_project in self.sample_projects.items():
            for name2, sample_id in sample_project.sample_ids.items():
                for sample in sample_id.samples:
                    for lane in sample.lane.split('+'):
                        lane2samples[lane].append(sample)
        for lane, samples in lane2samples.items():
            if len(samples) > 1:
                raise AnalysisDriverError('%s samples in lane %s despite has barcode set to %s'%(len(samples, lane, self.has_barcode)))

    def generate_mask(self, mask):
        """
        Translate:
            <Read IsIndexedRead=N Number=1 NumCycles=151/>
            <Read IsIndexedRead=Y Number=2 NumCycles=8/>
            <Read IsIndexedRead=N Number=3 NumCycles=151/>
        to 'y150n,i8,y150n'.
        """
        self.debug('Generating mask...')

        if self.has_barcode:
            barcode_len = self.check_barcodes()
            out = ['y' + str(mask.num_cycles(mask.upstream_read) - 1) + 'n']

            for i in mask.index_lengths:
                diff = i - barcode_len
                out.append('i' + str(barcode_len) + 'n' * diff)

            out.append('y' + str(mask.num_cycles(mask.downstream_read) - 1) + 'n')
            self.debug(out)
            return ','.join(out)

        elif not self.has_barcode:

            out = ['y' + str(mask.num_cycles(mask.upstream_read) - 1) + 'n']

            out.append('y' + str(mask.num_cycles(mask.downstream_read) - 1) + 'n')
            self.debug(out)
            return ','.join(out)

    def validate(self, mask):
        """
        Ensure that the SampleSheet is consistent with itself and RunInfo
        """
        self.debug('Validating...')
        if mask.has_barcodes:
           if self.check_barcodes() != mask.barcode_len:
            self.error(
                'Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' %
                (self.check_barcodes(), mask.barcode_len)
            )
            return False
        self.debug('Done. Now validating RunInfo')
        return True

    def get_samples(self, sample_project, sample_id):
        return self.sample_projects[sample_project].sample_ids[sample_id].samples

    def _populate(self):
        f, header = _read_sample_sheet(self.filename)
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

    def _get_sample_project(self, name):
        sample_project = ValueError('Could not add sample project ' + name)
        try:
            sample_project = self.sample_projects[name]
        except KeyError:
            sample_project = SampleProject(name)
            self.sample_projects[name] = sample_project
        finally:
            return sample_project

    @staticmethod
    def _get_column(header, name):
        possible_fields = config.sample_sheet_config['column_names'][name]
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


class SampleProject:
    """
    Represents a sample project, and contains a list of SampleID objects.
    """
    def __init__(self, name):
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
