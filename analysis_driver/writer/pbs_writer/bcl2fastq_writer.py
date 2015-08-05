__author__ = 'mwham'
import os.path
from analysis_driver.writer.pbs_writer import PBSWriter


class BCL2FastqWriter(PBSWriter):
    def __init__(self, pbs_name, job_name, log_file, walltime='24', cpus='12', mem='24', queue='uv2000'):
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file, queue)

    def _bcl2fastq(self, mask, input_dir, fastq_path):
        self.log('Writing bcl2fastq command')
        self.write_line(
            'bcl2fastq -l INFO --runfolder-dir %s --output-dir %s --sample-sheet %s --use-bases-mask %s\n' % (
                input_dir, fastq_path, os.path.join(input_dir, 'SampleSheet.csv'), mask
            )
        )

    def write(self, mask, input_dir, fastq_path):
        self._bcl2fastq(mask, input_dir, fastq_path)
        self.save()
