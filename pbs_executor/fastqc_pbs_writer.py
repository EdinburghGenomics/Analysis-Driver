from .pbs_writer import PBSWriter


class FastqcPBSWriter(PBSWriter):
    def __init__(self, pbs_name, job_name, log_file, walltime='6', cpus='8', mem='3'):
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file)

    def _fastqc(self, input_dir):
        self.write_line('FASTQ_FILES=`find ' + input_dir + ' -name \'*.fastq.gz\'`\n')
        self.write_line('fastqc --nogroup -t 8 -q $FASTQ_FILES\n')

    def write(self, input_dir, run_dir):
        self._fastqc(input_dir)
        self.write_line('touch ' + run_dir + '/.fastqc_complete')
        self.save()
