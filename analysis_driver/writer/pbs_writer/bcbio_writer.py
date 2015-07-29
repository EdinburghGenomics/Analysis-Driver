__author__ = 'mwham'
from .pbs_writer import PBSWriter
from analysis_driver import config


class BCBioWriter(PBSWriter):
    """
    Writes a PBS script to run BCBio alignment/variant calling
    """
    def __init__(self, pbs_name, job_name, log_file, walltime='72', cpus='8', mem='64', queue='uv2000'):
        """
        See superclass
        """
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file, queue)

    def _bcbio(self, bcbio_path, run_yaml, workdir, cores=16):
        """
        Write commands to run BCBio
        :param str bcbio_path: Path to bcbio_nextgen.py executable
        :param str run_yaml: Path to the BCBio project yaml
        :param str workdir: The directory to rrun BCBio in
        :param str cores: Number of cores for BCBio to use
        """
        self.info('Writing BCBio command')

        self.write_line('export JAVA_HOME=' + config.default['gatk']['java_home'])
        self.write_line('export JAVA_BINDIR=' + config.default['gatk']['java_bindir'])
        self.write_line('export JAVA_ROOT=' + config.default['gatk']['java_root'] + '\n')

        cmd = '%s %s -n %s --workdir %s\n' % (bcbio_path, run_yaml, cores, workdir)
        self.info(cmd)
        self.write_line(cmd)

    def write(self, bcbio_path, run_yaml, workdir):
        self._bcbio(bcbio_path, run_yaml, workdir)
        self.save()
