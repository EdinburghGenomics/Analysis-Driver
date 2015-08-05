__author__ = 'mwham'
from .pbs_writer import PBSWriter
from analysis_driver import config


class BCBioWriter(PBSWriter):
    """
    Writes a PBS script to run BCBio alignment/variant calling
    """
    def __init__(self, pbs_name, job_name, log_file, num_samples,
                 walltime='72', cpus='8', mem='64', queue='uv2000'):
        """
        See superclass
        """
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file, queue)
        self.num_samples = num_samples
        self._setup(self.num_samples)
        self.job_number = 0
    
    def write(self):
        self.write_line('*) echo "Unexpected PBS_ARRAY_INDEX: $PBS_ARRAY_INDEX"')
        self.write_line('esac')
        self.info('Job number: %s. Num samples: %s' % (self.job_number, self.num_samples))
        assert self.job_number == self.num_samples
        self.save()

    def _setup(self, num_samples, cores=16):
        """
        Write commands to run BCBio
        :param str bcbio_path: Path to bcbio_nextgen.py executable
        :param str run_yaml: Path to the BCBio project yaml
        :param str workdir: The directory to run BCBio in
        :param str cores: Number of cores for BCBio to use
        """
        self.info('Writing BCBio command')

        wt = self.write_line

        wt('#PBS -J 1-' + str(num_samples) + '\n')

        wt('export JAVA_HOME=' + config.default['gatk']['java_home'])
        wt('export JAVA_BINDIR=' + config.default['gatk']['java_bindir'])
        wt('export JAVA_ROOT=' + config.default['gatk']['java_root'] + '\n')
        
        wt('case $PBS_ARRAY_INDEX in')

    def add_bcbio_job(self, run_yaml, workdir, cores=16):
        self.job_number += 1
        cmd = '%s %s -n %s --workdir %s' % (config.default['bcbio'], run_yaml, cores, workdir)
        self.info(cmd)
        self.write_line(str(self.job_number) + ') ' + cmd)
        self.write_line(';;')

