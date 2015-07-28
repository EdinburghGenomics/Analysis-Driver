__author__ = 'mwham'
import csv
import os.path
from analysis_driver.util import AppLogger, fastq_handler


class BCBioCSVWriter(AppLogger):
    def __init__(self, fastq_dir, run_dir, sample_sheet):
        csv_file = os.path.join(run_dir, 'samples.csv')
        self.info('Csv file: ' + csv_file)
        self.samples = open(csv_file, 'w')

        self.writer = csv.writer(self.samples)
        self.sample_sheet = sample_sheet
        self.fastq_dir = fastq_dir

    def write(self):
        self.writer.writerow(['samplename', 'description', 'batch'])
        for name, sample_project in self.sample_sheet.sample_projects.items():
            fastqs = fastq_handler.find_fastqs(self.fastq_dir, name)
            for fq in fastqs:
                self.writer.writerow([fq, '', name])
        self.samples.close()
        # TODO: check that this works on Ultra
        self.info('Written csv file')

