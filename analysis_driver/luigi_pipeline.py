__author__ = 'mwham'
import luigi
import os.path
from glob import glob
from analysis_driver import executor, writer, reader, util, clarity
from analysis_driver.quality_control import genotype_validation
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('luigi')


class Stage(luigi.Task):

    input_run_folder = luigi.Parameter()
    job_dir = luigi.Parameter()

    lf_name = None

    @property
    def sample_sheet_path(self):
        return os.path.join(self.input_run_folder, 'SampleSheet_analysis_driver.csv')

    @property
    def run_id(self):
        return os.path.basename(self.input_run_folder)

    @property
    def job_dir(self):
        return os.path.join(cfg['jobs_dir'], self.run_id)

    @property
    def fastq_dir(self):
        return os.path.join(
            self.job_dir, 'fastq'
        )

    @staticmethod
    def touch(f):
        open(f, 'w').close()

    @property
    def lf(self):
        return os.path.join(self.job_dir, self.lf_name)

    def output(self):
        return luigi.LocalTarget(self.lf)


class Bcl2Fastq(Stage):
    lf_name = '.bcl2fastq.done'

    def run(self):
        sample_sheet = reader.SampleSheet(self.input_run_folder)
        run_info = reader.RunInfo(self.input_run_folder)
        mask = sample_sheet.generate_mask(run_info.mask)

        exit_status = executor.execute(
            [
                writer.bash_commands.bcl2fastq(
                    self.input_run_folder,
                    self.fastq_dir,
                    sample_sheet.filename,
                    mask
                )
            ],
            job_name='bcl2fastq',
            run_id=self.run_id,
            walltime=32,
            cpus=8,
            mem=32
        ).join()
        print(exit_status)
        self.touch(self.lf)
        return exit_status


class Fastqc(Stage):
    lf_name = '.fastqc.done'

    @staticmethod
    def requires():
        return Bcl2Fastq()

    def run(self):
        exit_status = executor.execute(
            [writer.bash_commands.fastqc(fq) for fq in util.fastq_handler.find_all_fastqs(self.fastq_dir)],
            job_name='fastqc',
            run_id=self.run_id,
            walltime=6,
            cpus=1,
            mem=2
        ).join()
        print(exit_status)
        self.touch(self.lf)
        return exit_status


'''class GenotypeValidation(Stage):

    lf_name = '.genotype_validation.done'

    @staticmethod
    def requires():
        return Bcl2Fastq()

    def run(self):
        fqs = {}
        for sample_project in data['sample_id_fastqs']:
            for sample_id in data['sample_id_fastqs'][sample_project]:
                fqs[sample_id] = data['sample_id_fastqs'][sample_project][sample_id]
        gen_val = genotype_validation.GenotypeValidation(fqs, self.run_id)
        gen_val.start()
        self.touch(self.lf)
        return gen_val.join()'''


class BCBio(Stage):
    sample_project = luigi.Parameter()
    sample_id = luigi.Parameter()
    fastqs = luigi.Parameter()

    lf_name = '.bcbio.done'

    @property
    def user_sample_id(self):
        u = clarity.get_user_sample_name(self.sample_id)
        if not u:
            u = self.sample_id
        return u

    @staticmethod
    def requires():
        return Bcl2Fastq()

    def prepare_samples(self):
        cmd = util.bcbio_prepare_samples_cmd(
            self.job_dir,
            self.sample_id,
            self.fastqs,
            user_sample_id=self.user_sample_id
        )
        return executor.execute([' '.join(cmd)], env='local').join()

    def prepare_template(self):
        run_template = os.path.join(
            os.path.dirname(__file__),
            '..', 'etc', 'bcbio_alignment_' + cfg['genome'] + '.yaml'
        )
        if not os.path.isfile(run_template):
            raise AnalysisDriverError(
                'Could not find BCBio run template \'%s\'. Is the correct genome set?' % run_template
            )

        original_dir = os.getcwd()
        os.chdir(self.job_dir)

        bcbio_dir = os.path.join(self.job_dir, 'samples_' + self.sample_id + '-merged')
        cmd = ' '.join(
            [
                os.path.join(cfg['bcbio'], 'bin', 'bcbio_nextgen.py'),
                '-w template',
                run_template,
                bcbio_dir,
                bcbio_dir + '.csv'
            ] + self.fastqs
        )

        exit_status = executor.execute([cmd], env='local').join()
        os.chdir(original_dir)
        return exit_status

    def run_bcbio(self):

        cmd = writer.bash_commands.bcbio(
            os.path.join(
                self.job_dir,
                'samples_' + self.sample_id + '-merged',
                'config',
                'samples_' + self.sample_id + '-merged.yaml'
            ),
            os.path.join(
                self.job_dir,
                'samples_' + self.sample_id + '-merged',
                'work'
            ),
            threads=10
        )
        return executor.execute(
            [cmd],
            prelim_cmds=writer.bash_commands.bcbio_env_vars(),
            job_name='bcbio',
            run_id=self.run_id,
            walltime=96,
            cpus=10,
            mem=64
        ).join()

    def run(self):
        exit_status = 0
        for stage in (self.prepare_samples, self.prepare_template, self.run_bcbio):
            ex = stage()
            exit_status += ex
            if exit_status:
                break
        self.touch(self.lf)
        return exit_status


class DivergeSamples(Stage):

    lf_name = '.diverge_samples.done'

    @staticmethod
    def requires():
        for sample_project in data['sample_projects']:
            for sample_id in data['sample_projects'][sample_project]:
                yield BCBio(sample_project=sample_project, sample_id=sample_id)

    def output(self):
        return self.input()

    def run(self):
        self.touch(self.lf)


class DataOutput(Stage):

    @staticmethod
    def requires():
        return DivergeSamples()

    def run(self):
        for req in self.requires().requires():
            yield req

            exit_status = self.execute(req)
            app_logger.info('Finished transfer of %s with exit status %s' % (req.sample_id, exit_status))

    def execute(self, req):

        output_dir = os.path.join(cfg['output_dir'], req.sample_project, req.sample_id)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        bcbio_source_dir = os.path.join(
            self.job_dir,
            'samples_' + req.sample_id + '-merged',
            'final',
            req.sample_id
        )
        merged_fastq_dir = os.path.join(
            self.job_dir,
            'merged'
        )
        source_path_mapping = {
            'vcf': bcbio_source_dir,
            'bam': bcbio_source_dir,
            'fastq': merged_fastq_dir
        }
        app_logger.info('Beginning output transfer')
        return util.transfer_output_files(req.sample_id, output_dir, source_path_mapping)


def pipeline(input_run_folder):
    reader.transform_sample_sheet(input_run_folder)
    sample_sheet = reader.SampleSheet(input_run_folder)
    run_info = reader.RunInfo(input_run_folder)

    data['mask'] = sample_sheet.generate_mask(run_info.mask)

    data['sample_id_fastqs'] = {}
    for sample_project in sample_sheet.sample_projects:
        data['sample_id_fastqs'][sample_project] = {}
        for sample_id in sample_sheet.sample_projects[sample_project].sample_ids:

            data['sample_projects'][sample_project].append(sample_id)

    print(data)

    luigi.run(
        cmdline_args=[
            '--input-run-folder', input_run_folder,
            '--Stage-input-run-folder', input_run_folder
        ],
        main_task_cls=DataOutput,
        local_scheduler=True
    )
