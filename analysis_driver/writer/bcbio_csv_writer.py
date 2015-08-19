__author__ = 'mwham'
import csv
import os.path
from logging import getLogger
from analysis_driver.util import AppLogger


class BCBioCSVWriter(AppLogger):
    """
    Writes a BCBio csv sample file based on an input SampleSheet object.
    """
    def __init__(self, run_dir, sample_id, fastqs):
        """
        :param str sample_id: A sample id to assign to the csv
        :param str run_dir: Full path to a run folder
        :param list fastqs: Fastqs to assign to the sample id
        """
        self.csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
        self.sample_id = sample_id
        self.fastqs = fastqs
        self.info('csv file: ' + self.csv_file)
        self.samples = open(self.csv_file, 'w')
        self.writer = csv.writer(self.samples)

    def write(self):
        """
        Iterate through sample_projects in self.sample_sheet. For each, find relevant fastq files and write a
        line to self.samples mapping the fastq name to the sample project id.
        """
        self.writer.writerow(['samplename', 'description'])
        sample_id = self.sample_id
        for fq in self.fastqs:
            self.writer.writerow([fq, sample_id])
        self.samples.close()
        self.info('Written csv file')


logger = getLogger(__name__)


def write_bcbio_csv(run_dir, sample_id, fastqs):
    csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)

        writer.writerow(['samplename', 'description'])
        for fq in fastqs:
            writer.writerow([fq, sample_id])

    return csv_file
