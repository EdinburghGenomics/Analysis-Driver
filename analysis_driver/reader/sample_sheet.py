__author__ = 'mwham'
import csv
import os
from analysis_driver.util import AppLogger


class SampleSheet(AppLogger):
    """
    Represents an instance of SampleSheet.csv. It ignores all lines until '[Data]' is found in the first
    column, then starts reading the CSV.
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
        self._populate()
        self.debug('Sample project entries: ' + str(self.sample_projects))

    def _populate(self):
        reader = csv.DictReader(self.file)
        for row in reader:
            if any(row):
                sample_project = row['Sample_Project']
                sample_id = row['Sample_ID']

                new_sample = Sample(
                    sample_project=sample_project,
                    lane=row['Lane'],
                    id=sample_id,
                    name=row['Sample_Name'],
                    index=row['Index']  # ,
                    # index2=row['Index2'],
                    # plate=row['Sample_Plate'],
                    # well=row['Sample_Well'],
                    # genome_folder=row['GenomeFolder']
                )

                sample_project_obj = self._get_sample_project(sample_project)
                sample_id_obj = sample_project_obj.get_sample_id(sample_id)
                sample_id_obj.add_sample(new_sample)

    def check_barcodes(self):
        """
        For each sample project, check that all the DNA barcodes are the same length
        :return: The DNA barcode length for the sample sheet
        :rtype: int
        """
        last_sample = None
        for name, sample_project in self.sample_projects.items():
            for name2, sample_id in sample_project.sample_ids.items():
                for sample in sample_id.samples:
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
        self.fastq_files = []

    def add_sample(self, sample):
        assert sample.sample_project == self.sample_project and sample.id == self.name,\
            'Adding invalid sample project to ' + self.name + ': ' + sample.sample_project
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
