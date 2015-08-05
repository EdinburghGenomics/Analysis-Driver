__author__ = 'mwham'
import csv
import os.path
from analysis_driver.util import AppLogger, fastq_handler


class BCBioSamplePrep(AppLogger):
    def __init__(self, fastq_dir, run_dir, sample_sheet):
        """
        :param str fastq_dir: Full path to a dir to search for fastq files
        :param str run_dir: Full path to a run folder
        :param SampleSheet sample_sheet: A SampleSheet object containing data on samples to assign
        """
        self.fastq_dir = fastq_dir
        self.run_dir = run_dir
        self.sample_sheet = sample_sheet

    def write(self):
        csvs = []

        for name, sample_project in self.sample_sheet.sample_projects.items():
            fastqs = fastq_handler.find_fastqs(self.fastq_dir, name)

            for sample_id in sample_project.sample_ids:
                self.info('Writing csv file for: ' + sample_id)
                writer = BCBioCSVWriter(self.run_dir, sample_id, fastqs[sample_id])
                writer.write()
                csvs.append((sample_id, writer.csv_file, fastqs[sample_id]))
        return csvs


class BCBioCSVWriter(AppLogger):
    """
    Writes a BCBio csv sample file based on an input SampleSheet object.
    """
    def __init__(self, run_dir, sample_id, fastqs):
        """
        :param str fastq_dir: Full path to a dir to search for fastq files
        :param str run_dir: Full path to a run folder
        :param SampleSheet sample_sheet: A SampleSheet object containing data on samples to assign
        """
        self.csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
        self.sample_id = sample_id
        self.fastqs = fastqs
        self.info('Csv file: ' + self.csv_file)
        self.samples = open(self.csv_file, 'w')
        self.writer = csv.writer(self.samples)

    def write(self):
        """
        Iterate through sample_projects in self.sample_sheet. For each, find relevant fastq files and write a
        line to self.samples mapping the fastq name to the sample project id.
        """
        self.writer.writerow(['samplename', 'description'])
        id = self.sample_id
        for fq in self.fastqs:
            self.writer.writerow([fq, id])
        self.samples.close()
        self.info('Written csv file')

