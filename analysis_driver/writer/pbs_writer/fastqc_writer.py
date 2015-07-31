__author__ = 'mwham'
from analysis_driver.writer.pbs_writer import PBSWriter


class FastqcWriter(PBSWriter):
    """
    Writes a PBS script to run fastqc
    """
    def __init__(self, pbs_name, job_name, log_file, fastqs=None, walltime='6', cpus='8', mem='3', queue='uv2000'):
        """
        See superclass
        """
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file, queue)
        self.fastqs = fastqs

    def _fastqc(self, fastqs):
        """
        Write commands to run fastqc
        :param input_dir: Path to a directory containing input .fastq files
        """
        self.info('Writing fastqc command')
        wt = self.write_line

        wt('#PBS -J 1-' + str(len(fastqs)) + '\n')

        wt('case $PBS_ARRAY_INDEX in')
        for idx, fastq in enumerate(fastqs):
            wt('%s) fastqc --nogroup -q %s' % (self._shidx(idx), fastq))
            wt(';;')
        wt('*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"')
        wt('esac')

    def write(self):
        """
        Create a .fastqc_complete lock file to tell the AnalysisDriver when to proceed to later steps
        """
        self._fastqc(self.fastqs)
        self.save()

    @staticmethod
    def _shidx(idx):
        return str(idx + 1)
