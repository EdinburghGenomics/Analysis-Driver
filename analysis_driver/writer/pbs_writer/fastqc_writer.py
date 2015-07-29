__author__ = 'mwham'
from analysis_driver.writer.pbs_writer import PBSWriter


class FastqcWriter(PBSWriter):
    """
    Writes a PBS script to run fastqc
    """
    def __init__(self, pbs_name, job_name, log_file, walltime='6', cpus='8', mem='3', queue='uv2000'):
        """
        See superclass
        """
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file, queue)

    def _fastqc(self, input_dir):
        """
        Write commands to run fastqc
        :param input_dir: Path to a directory containing input .fastq files
        """
        # TODO: find fastqs through Python rather than Bash
        self.info('Writing fastqc command')
        self.write_line('FASTQ_FILES=`find ' + input_dir + ' -name \'*.fastq.gz\'`\n')
        self.write_line('fastqc --nogroup -t 8 -q $FASTQ_FILES\n')

    def write(self, input_dir, run_dir):
        """
        Create a .fastqc_complete lock file to tell the AnalysisDriver when to proceed to later steps
        """
        self._fastqc(input_dir)
        self.write_line('touch ' + run_dir + '/.fastqc_complete')
        self.save()
