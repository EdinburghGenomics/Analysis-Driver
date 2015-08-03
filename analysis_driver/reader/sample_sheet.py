__author__ = 'mwham'
import csv
import os
from analysis_driver.util import AppLogger


class SampleSheet(AppLogger):
    """
    Represents an instance of SampleSheet.csv. It ignores all lines until '[Data]' is found in the first
    column, then starts reading the CSV.

    Public properties:
        sample_projects: A dict mapping sample project names to SampleProject objects (dict[str, SampleProject])

    Public methods:
        check_barcodes()
    """
    def __init__(self, data_dir):
        """
        :param str data_dir:  A file path to the input_data folder containing SampleSheet.csv
        """
        self.file = open(os.path.join(data_dir, 'SampleSheet.csv'), 'r')
        # read lines until [Data] marker
        counter = 1
        while not next(self.file).startswith('[Data]'):
            pass  # Do nothing until [Data]
            counter += 1
        self.debug('Starting reading sample sheet at line ' + str(counter))

        self.sample_projects = {}  # {name: samples} {str: Sample}
        self._populate(self.sample_projects)
        self.debug('Sample project entries: ' + str(self.sample_projects))

    def _populate(self, samples):
        reader = csv.DictReader(self.file)
        for row in reader:
            if any(row):
                try:
                    sample_project = row['Project_Name']
                except KeyError:
                    sample_project = row['Sample_Project']
                new_sample = Sample(
                    sample_project=sample_project,
                    lane=row['Lane'],
                    id=row['Sample_ID'],
                    name=row['Sample_Name'],
                    index=row['Index'],
                    index2=row['Index2'],
                    plate=row['Sample_Plate'],
                    well=row['Sample_Well'],
                    genome_folder=row['GenomeFolder']
                )
                try:
                    samples[sample_project].add_sample(new_sample)
                except KeyError:
                    samples[sample_project] = SampleProject(sample_project, new_sample)

    def check_barcodes(self):
        """
        For each sample project, check that all the DNA barcodes are the same length
        :return: The DNA barcode length for the sample sheet
        :rtype: int
        """
        last_sample = None
        for name, sample_project in self.sample_projects.items():
            for sample in sample_project.samples:
                try:
                    if len(sample.index) == len(last_sample.index):
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

        return len(last_sample.index)


class SampleProject:
    """
    Represents a sample project, and contains a list of Sample objects.

    Public properties:
        name: The sample project name
        samples: A list of Sample objects

    Public methods:
        add_sample()
    """
    def __init__(self, name, new_sample):
        """
        :param name: The saple project name
        :param new_sample: The first Sample to initialise the SampleProject with
        """
        self.name = name
        self.samples = [new_sample]

    def add_sample(self, sample):
        """
        Append a Sample object to self.samples
        :param sample: A Sample object to append
        """
        if sample.sample_project != self.name:
            raise AssertionError(
                'Adding invalid sample project to ' + self.name + ': ' + sample.sample_project
            )
        else:
            self.samples.append(sample)


class Sample:
    """
    This represents a Sample, i.e. a line in SampleSheet.csv below the '[Data]' marker. Supports dict-style
    attribute fetching, e.g. sample_object['lane']
    """
    def __init__(self, **kwargs):
        """
        :param kwargs: The columns of SampleSheet.csv and corresponding values to assign
        :raises: AssertionError if  a column is invalid
        """
        self.data = {}
        for k, v in kwargs.items():
            assert k in [
                'sample_project', 'lane', 'id', 'name', 'index', 'index2', 'plate', 'well',
                'genome_folder'
            ]
            self.data[k] = v

    def __getattr__(self, attr):
        return self.data[attr]
