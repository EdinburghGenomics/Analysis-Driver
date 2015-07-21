import os.path
from .pbs_writer import PBSWriter
import util


class BCBioPBSWriter(PBSWriter):
    def __init__(self, pbs_name, job_name, log_file, walltime='72', cpus='8', mem='64'):
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file)

    @staticmethod
    def get_fastqs(fastq_dir, sample_project):
        f = []
        fastqs = os.path.join(fastq_dir, sample_project)
        for sample_id in os.listdir(fastqs):
            f = f + [
                os.path.join(
                    fastqs, sample_id, fq
                ) for fq in os.listdir(
                    os.path.join(fastqs, sample_id)
                ) if fq.endswith('.fastq.gz')
            ]
        print('Fastqs: ', f)
        return f

    @staticmethod
    def setup_bcbio_run(bcbio, template, csv_file, run_dir, fastqs):
        util.localexecute(
            bcbio,
            '-w',
            'template',
            template,
            run_dir,
            csv_file,
            *fastqs
        )

    def _bcbio(self, bcbio_path, run_yaml, workdir, cores=16):
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

