__author__ = 'mwham'
import os.path
from analysis_driver.writer.pbs_writer import PBSWriter


class BCL2FastqWriter(PBSWriter):
    """
    Writes a PBS script to run bcl2fastq
    """
    def __init__(self, pbs_name, job_name, log_file, walltime='32', cpus='12', mem='24', queue='uv2000'):
        """
        See superclass
        """
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file, queue)

    def _bcl2fastq(self, mask, input_dir, fastq_path):
        """
        Write commands to run bcl2fastq
        :param str mask: A mask to use, as generated by reader.RunInfo
        :param input_dir: Path to the input dir containing the bcl files
        :param fastq_path: Path to the dir in which to send generated .fastqs
        """
        self.info('Writing bcl2fastq command')
        self.write_line(
            'bcl2fastq -l INFO --runfolder-dir %s --output-dir %s --sample-sheet %s --use-bases-mask %s\n' % (
                input_dir, fastq_path, os.path.join(input_dir, 'SampleSheet.csv'), mask
            )
        )

    def write(self, job_dir, mask, input_dir, fastq_path):
        self._bcl2fastq(mask, input_dir, fastq_path)
        self.write_line('touch ' + job_dir + '/.bcl2fastq_complete')

        self.save()
