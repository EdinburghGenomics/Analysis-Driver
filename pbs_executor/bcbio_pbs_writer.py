import os.path
from .pbs_writer import PBSWriter
import util


class BCBioPBSWriter(PBSWriter):
    def __init__(self, pbs_name, job_name, log_file, walltime='72', cpus='8', mem='64'):
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file)

    @staticmethod
    def get_fastqs(fastq_dir, sample_project):
        print('[BCBioPBSWriter]', fastq_dir, sample_project)
        fastqs = os.path.join(fastq_dir, sample_project)
        return [os.path.join(fastqs, f) for f in os.listdir(fastqs) if f.endswith('fastq.gz')]

    @staticmethod
    def setup_bcbio_run(bcbio, csv_file, template, sample_project, fastqs):
        util.localexecute(bcbio, '-w', 'template', csv_file, template, sample_project, *fastqs)

    def _bcbio(self, bcbio_path, run_yaml, workdir, cores=16):
        # Java paths. TODO: get these into the YAML config
        self.write_line('export JAVA_HOME=/home/U008/edingen/Applications/jdk1.7.0_76/')
        self.write_line('export JAVA_BINDIR=/home/U008/edingen/Applications/jdk1.7.0_76/bin')
        self.write_line('export JAVA_ROOT=/home/edingen/Applications/jdk1.7.0_76/\n')

        self.write_line('%s %s -n %s --workdir %s\n' % (bcbio_path, run_yaml, cores, workdir))

    def write(self, bcbio_path, run_yaml, workdir):
        self._bcbio(bcbio_path, run_yaml, workdir)
        self.save()
