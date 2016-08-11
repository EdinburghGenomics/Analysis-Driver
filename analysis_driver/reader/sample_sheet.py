import csv
from os.path import join
from collections import defaultdict
from egcg_core.app_logging import AppLogger, logging_default as log_cfg
from egcg_core.clarity import get_library_id
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import sample_sheet_config

app_logger = log_cfg.get_logger('reader')


def transform_sample_sheet(data_dir, remove_barcode=False):
    """
    Read SampleSheet.csv, translate column names and write to SampleSheet_analysis_driver.csv
    :param str data_dir: Full path to a data directory containing SampleSheet.csv
    :param bool remove_barcode: Whether to remove barcodes from the sample sheet
    """
    original, header = _read_sample_sheet(join(data_dir, 'SampleSheet.csv'))
    transformations = sample_sheet_config.get('transformations', [])
    old_col_names = original.readline().strip().split(',')
    new_col_names = [transformations.get(x, x) for x in old_col_names]

    out_lines = [line.strip() for line in header] + ['[Data],', ','.join(new_col_names)]

    if remove_barcode:
        for line in original:
            sp_line = line.strip().split(',')
            sp_line[new_col_names.index('Index')] = ''
            sp_line[new_col_names.index('Index2')] = ''
            out_lines.append(','.join(sp_line))
    else:
        out_lines.extend(line.strip() for line in original)

    original.close()  # close the open file from _read_sample_sheet
    out_lines.append('')
    with open(join(data_dir, 'SampleSheet_analysis_driver.csv'), 'w') as new_sample_sheet:
        new_sample_sheet.write('\n'.join(out_lines))


def _read_sample_sheet(sample_sheet):
    """
    Scan down a sample sheet until a [Data] line, then return the (open!) file object for further reading.
    :return tuple[file, list] f, header: The sample sheet file object, and all lines above [Data]
    """
    f = open(sample_sheet, 'r')
    app_logger.debug('Opened ' + f.name)
    counter = 1
    header = []
    for line in f:
        if line.startswith('[Data]'):
            app_logger.debug('Reading sample sheet from line ' + str(counter))
            return f, header
        else:
            counter += 1
            header.append(line)
    f.close()
    return None, None


class SampleSheet(AppLogger):
    def __init__(self, filename, has_barcode=True):
        self.projects = {}  # {name: samples} {str: Sample}
        self.filename = filename
        self._populate()
        self.debug('Projects: ' + str(self.projects))
        self.has_barcode = has_barcode
        if not self.has_barcode:
            self._validate_one_sample_per_lane()

    def check_barcodes(self):
        """
        For each project in the sample sheet, check that all the DNA barcodes are the same length
        :return: The DNA barcode length for the sample sheet
        """
        last_barcode_len = None
        for project_id, project_obj in self.projects.items():
            self.debug('Checking project ' + project_id)
            for sample_id, sample_obj in project_obj.sample_ids.items():
                self.debug('Checking sample id ' + sample_id)
                for sample in sample_obj.samples:
                    if last_barcode_len and last_barcode_len != len(sample.barcode):
                        raise AnalysisDriverError(
                            "Odd barcode length for %s: %s (%s) in project '%s' " % (
                                sample.id, sample.barcode, len(sample.barcode), project_id
                            )
                        )
                    last_barcode_len = len(sample.barcode)
        self.debug('Barcode check done. Barcode len: %s', last_barcode_len)
        return last_barcode_len

    def _validate_one_sample_per_lane(self):
        """
        Check that only one sample is present in each lane and raise AnalysisDriverError if more than one is
        found.
        """
        samples_per_lane = defaultdict(list)
        for p in self.projects.values():
            for sample_id in p.sample_ids.values():
                for sample in sample_id.samples:
                    for lane in sample.lane.split('+'):
                        samples_per_lane[lane].append(sample)

        for lane, samples in samples_per_lane.items():
            if len(samples) > 1:
                raise AnalysisDriverError(
                    'Non-barcoded run, but lane %s has %s samples.' % (lane, len(samples))
                )

    def generate_mask(self, mask):
        """
        Translate:
            <Read IsIndexedRead=N Number=1 NumCycles=151/>
            <Read IsIndexedRead=Y Number=2 NumCycles=8/>
            <Read IsIndexedRead=N Number=3 NumCycles=151/>
        to 'y150n,i8,y150n', depending on self.has_barcode.
        :param .run_info.Mask mask: A Mask object extracted from RunInfo.xml
        """
        self.debug('Generating mask')
        out = ['y' + str(mask.num_cycles(mask.upstream_read) - 1) + 'n']

        if self.has_barcode:
            barcode_len = self.check_barcodes()
            for i in mask.index_lengths:
                diff = i - barcode_len
                out.append('i' + str(barcode_len) + 'n' * diff)

        out.append('y' + str(mask.num_cycles(mask.downstream_read) - 1) + 'n')
        self.debug(out)
        return ','.join(out)

    def validate(self, mask):
        """
        Ensure that the SampleSheet is consistent with itself and RunInfo.
        :param .run_info.Mask mask: Mask object to check against
        """
        self.debug('Validating...')
        if mask.has_barcodes and self.check_barcodes() != mask.barcode_len:
            self.error(
                'Barcode mismatch: %s (SampleSheet.csv) and %s (RunInfo.xml)' %
                (self.check_barcodes(), mask.barcode_len)
            )
            return False
        return True

    def get_samples(self, project_id, sample_id):
        return self.projects[project_id].sample_ids[sample_id].samples

    def _populate(self):
        f, header = _read_sample_sheet(self.filename)
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        counter = 0
        for line in reader:
            if any(line):
                counter += 1
                project_id = line[self._get_column(cols, 'project_id')]
                sample_id = line[self._get_column(cols, 'sample_id')]

                new_sample = Sample(
                    project_id=project_id,
                    lane=line[self._get_column(cols, 'lane')],
                    sample_id=sample_id,
                    default_library_id=line[self._get_column(cols, 'library_id')],
                    barcode=line[self._get_column(cols, 'barcode')],
                    **self._get_all_cols(
                        line,
                        ignore=('project_id', 'sample_id', 'library_id' 'lane', 'barcode')
                    )
                )

                project_obj = self._get_project(project_id)
                sample_id_obj = project_obj.get_sample_id(sample_id)
                sample_id_obj.add_sample(new_sample)
        f.close()
        self.debug('Added %s samples', counter)

    def _get_project(self, name):
        if name not in self.projects:
            self.projects[name] = Project(name)
        return self.projects[name]

    @staticmethod
    def _get_column(header, name):
        for f in sample_sheet_config['column_names'][name]:
            if f in header:
                return f

    @staticmethod
    def _get_all_cols(line_dict, ignore=None):
        return dict((k, v) for k, v in line_dict.items() if k not in ignore)


class Project:
    """Represents a project in the sample sheet, containing a list of SampleID objects"""
    def __init__(self, name):
        self.name = name
        self.sample_ids = {}

    def get_sample_id(self, name):
        if name not in self.sample_ids:
            self.sample_ids[name] = SampleID(name, project_id=self.name)
        return self.sample_ids[name]


class SampleID:
    def __init__(self, name, project_id):
        self.name = name
        self.project_id = project_id
        self.samples = []

    def add_sample(self, sample):
        assert sample.project_id == self.project_id and sample.sample_id == self.name,\
            'Adding invalid project id to ' + self.name + ': ' + sample.project_id
        self.samples.append(sample)


class Sample:
    """Represents a line in SampleSheet.csv below the '[Data]' marker."""
    _library_id = None

    def __init__(self, project_id, sample_id, default_library_id, lane, barcode, **kwargs):
        self.project_id = project_id
        self.sample_id = sample_id
        self.default_library_id = default_library_id
        self.lane = lane
        self.barcode = barcode
        self.extra_data = kwargs

    @property
    def library_id(self):
        """Query the Lims for the sample's library ID. Default to what is in the sample sheet."""
        if self._library_id is None:
            self._library_id = get_library_id(self.sample_id) or self.default_library_id
        return self._library_id
