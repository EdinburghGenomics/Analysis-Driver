__author__ = 'mwham'
import csv
import os.path


class BCBioCSVWriter:
    def __init__(self, fastq_dir, run_dir, sample_sheet):
        self.writer = csv.writer(open(os.path.join(run_dir, 'bcbio_samples.csv')))
        self.sample_sheet = sample_sheet
        self.fastq_dir = fastq_dir

    def write(self):
        self.writer.writerow(['samplename', 'description', 'batch', 'phenotype', 'sex', 'variant_regions'])
        for name, sample_project in self.sample_sheet.sample_projects.items():
            fastqs = self._find_fastqs(self.fastq_dir, name)
            for fq in fastqs:
                self.writer.writerow([fq, '', name, '', '', ''])

    @staticmethod
    def _find_fastqs(location, sample_project):
        fastqs = []
        fastq_dir = os.path.join(location, sample_project)
        sample_ids = [os.path.join(fastq_dir, x) for x in os.listdir(fastq_dir)]
        for id in sample_ids:
            fastqs = fastqs + [os.path.join(id, fq) for fq in os.listdir(id)]

        return fastqs
