__author__ = 'mwham'
from .pbs_writer import PBSWriter


class BCBioPBSWriter(PBSWriter):
    def __init__(self, pbs_name, job_name, log_file, walltime='72', cpus='8', mem='64'):
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file)

    def _bcbio(self, bcbio_path, run_yaml, workdir, cores=16):
        self.log('Writing BCBio command')
        # Java paths. TODO: get these into the YAML config
        self.write_line('export JAVA_HOME=/home/U008/edingen/Applications/jdk1.7.0_76/')
        self.write_line('export JAVA_BINDIR=/home/U008/edingen/Applications/jdk1.7.0_76/bin')
        self.write_line('export JAVA_ROOT=/home/edingen/Applications/jdk1.7.0_76/\n')

        cmd = '%s %s -n %s --workdir %s\n' % (bcbio_path, run_yaml, cores, workdir)
        self.log(cmd)
        self.write_line(cmd)

    def write(self, bcbio_path, run_yaml, workdir):
        self._bcbio(bcbio_path, run_yaml, workdir)
        self.save()
