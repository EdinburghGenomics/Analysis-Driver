__author__ = 'mwham'
import csv
import os.path
from analysis_driver.util import AppLogger, fastq_handler


class BCBioCSVWriter(AppLogger):
    """
    Writes a BCBio csv sample file based on an input SampleSheet object.
    """
    def __init__(self, fastq_dir, run_dir, sample_sheet):
        """
        :param str fastq_dir: Full path to a dir to search for fastq files
        :param str run_dir: Full path to a run folder
        :param SampleSheet sample_sheet: A SampleSheet object containing data on samples to assign
        """
        csv_file = os.path.join(run_dir, 'samples.csv')
        self.info('Csv file: ' + csv_file)
        self.samples = open(csv_file, 'w')

        self.writer = csv.writer(self.samples)
        self.sample_sheet = sample_sheet
        self.fastq_dir = fastq_dir

    def write(self):
        """
        Iterate through sample_projects in self.sample_sheet. For each, find relevant fastq files and write a
        line to self.samples mapping the fastq name to the sample project id.
        """
        self.writer.writerow(['samplename', 'description'])
        for name, sample_project in self.sample_sheet.sample_projects.items():
            fastqs = fastq_handler.find_fastqs(self.fastq_dir, name)
            print(fastqs)

            for sample_id in sample_project.sample_ids:
                for fq in fastqs[sample_id]:
                    self.writer.writerow([fq, sample_id])
        self.samples.close()
        self.info('Written csv file')

# TODO: write one csv per sample, pass through bcbio_prepare_samples
