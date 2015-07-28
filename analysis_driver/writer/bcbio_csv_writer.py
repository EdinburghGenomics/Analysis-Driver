__author__ = 'mwham'
import csv
import os.path
from analysis_driver.util import AppLogger


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
            fastqs = self._find_fastqs(self.fastq_dir, name)
            for fq in fastqs:
                self.writer.writerow([fq, '', name])
        self.samples.close()
        self.info('Written csv file')

    def _find_fastqs(self, location, sample_project):
        # TODO: push this up to util
        # TODO: have this search per sample within a sample project
        fastqs = []
        fastq_dir = os.path.join(location, sample_project)

        sample_ids = [os.path.join(fastq_dir, x) for x in self._list_dirs(fastq_dir)]
        for sample_id in sample_ids:
            fastqs = fastqs + [
                os.path.join(sample_id, fq) for fq in os.listdir(sample_id) if fq.endswith('.fastq.gz')
            ]

        self.info('Found ' + str(len(fastqs)) + 'fastq files')
        return fastqs

    @staticmethod
    def _list_dirs(d):
        return [x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x))]
