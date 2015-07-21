__author__ = 'mwham'
import csv
import os

from analysis_driver.util import AppLogger


class SampleSheet(AppLogger):
    def __init__(self, data_dir):
        self.file = open(os.path.join(data_dir, 'SampleSheet.csv'), 'r')
        # read lines until [Data] marker
        counter = 1
        while not next(self.file).startswith('[Data]'):
            pass  # Do nothing until [Data]
            counter += 1
        self.log('Start reading sample sheet at line ' + str(counter), 'DEBUG')

        self.sample_projects = {}  # {name: samples} {str: Sample}
        self._populate(self.sample_projects)
        self.log('Sample project entries: ' + str(self.sample_projects), 'DEBUG')

    def _populate(self, samples):
        reader = csv.DictReader(self.file)
        for row in reader:
            if any(row):
                sample_project = row['Sample_Project']
                new_sample = Sample(
                    sample_project=sample_project,
                    lane=row['Lane'],
                    id=row['Sample_ID'],
                    name=row['Sample_Name'],
                    index_id=row['I7_Index_ID'],
                    barcode=row['index'],
                    plate=row['Sample_Plate'],
                    well=row['Sample_Well'],
                    description=row['Description']
                )
                try:
                    samples[sample_project].add_sample(new_sample)
                except KeyError:
                    samples[sample_project] = SampleProject(sample_project, new_sample)

    def check_barcodes(self):
        """
        For each sample project, check that all the DNA barcodes are the same length
        TODO: yield the length for each sample project
        :return: The DNA barcode length for each sample project
        :rtype: int
        """
        for name, sample_project in self.sample_projects.items():
            last_sample = None
            for sample in sample_project.samples:
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

        return len(last_sample.barcode)


class SampleProject:
    def __init__(self, name, new_sample):
        self.name = name
        self.samples = [new_sample]

    def add_sample(self, sample):
        if sample.sample_project != self.name:
            raise AssertionError
        else:
            self.samples.append(sample)


class Sample:
    def __init__(self, **kwargs):
        self.data = {}
        for k, v in kwargs.items():
            assert k in [
                'sample_project', 'lane', 'id', 'name', 'index_id', 'barcode', 'plate', 'well',
                'description'
            ]
            self.data[k] = v

    def __getattr__(self, attr):
        return self.data[attr]
