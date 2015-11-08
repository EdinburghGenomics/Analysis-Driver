__author__ = 'mwham'
import luigi
import yaml
import os.path
from glob import glob
from analysis_driver import executor, writer, reader, util, clarity
# from analysis_driver.quality_control import genotype_validation
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('luigi')


# noinspection PyTypeChecker
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

    def execute(self):
        raise NotImplementedError

    def run(self):
        exit_status = self.execute()
        self.touch(self.lf)
        return exit_status


# noinspection PyTypeChecker
class Bcl2Fastq(Stage):
    lf_name = '.bcl2fastq.done'

    def execute(self):
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
        return exit_status


class Fastqc(Stage):
    lf_name = '.fastqc.done'

    def requires(self):
        return Bcl2Fastq()

    def execute(self):
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

    def execute(self):
        fqs = {}
        for sample_project in data['sample_id_fastqs']:
            for sample_id in data['sample_id_fastqs'][sample_project]:
                fqs[sample_id] = data['sample_id_fastqs'][sample_project][sample_id]
        gen_val = genotype_validation.GenotypeValidation(fqs, self.run_id)
        gen_val.start()
        self.touch(self.lf)
        return gen_val.join()'''


# noinspection PyTypeChecker
class BCBio(Stage):
    sample_project = luigi.Parameter()
    sample_id = luigi.Parameter()

    @property
    def lf_name(self):
        return '.bcbio_' + self.sample_id + '.done'

    @property
    def user_sample_id(self):
        u = clarity.get_user_sample_name(self.sample_id)
        if not u:
            u = self.sample_id
        return u

    @property
    def merged_dir(self):
        return os.path.join(self.job_dir, 'samples_' + self.sample_id + '-merged')

    @property
    def run_yaml(self):
        return os.path.join(self.merged_dir, 'config', 'samples_' + self.sample_id + '-merged.yaml')

    @property
    def fastqs(self):
        return glob(os.path.join(self.merged_dir, self.user_sample_id + '_R?.fastq.gz'))

    def requires(self):
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

        if exit_status:
            return exit_status

        with open(self.run_yaml, 'r') as i:
            run_config = yaml.load(i)
        run_config['fc_name'] = self.user_sample_id
        with open(self.run_yaml, 'w') as o:
            o.write(yaml.safe_dump(run_config, default_flow_style=False))

        return exit_status

    def run_bcbio(self):
        original_dir = os.getcwd()
        os.chdir(self.job_dir)

        cmd = writer.bash_commands.bcbio(self.run_yaml, workdir=os.path.join(self.merged_dir, 'work'), threads=10)
        exit_status = executor.execute(
            [cmd],
            prelim_cmds=writer.bash_commands.bcbio_env_vars(),
            job_name='bcbio',
            run_id=self.run_id,
            walltime=96,
            cpus=10,
            mem=64
        ).join()

        os.chdir(original_dir)
        return exit_status

    def execute(self):
        exit_status = 0
        for stage in (self.prepare_samples, self.prepare_template, self.run_bcbio):
            ex = stage()
            exit_status += ex
            if exit_status:
                break
        self.touch(self.lf)
        return exit_status


# noinspection PyTypeChecker
class DataOutput(Stage):

    @property
    def lf_name(self):
        return '.data_output_' + self.sample_id + '.done'

    sample_project = luigi.Parameter()
    sample_id = luigi.Parameter()

    @property
    def user_sample_id(self):
        u = clarity.get_user_sample_name(self.sample_id)
        if not u:
            u = self.sample_id
        return u

    def requires(self):
        return BCBio(sample_project=self.sample_project, sample_id=self.sample_id)

    def execute(self):
        exit_status = 0

        output_loc = os.path.join(self.output_dir, self.sample_project, self.sample_id)
        if not os.path.isdir(output_loc):
            os.makedirs(output_loc)

        for output_record in cfg['output_files']:
            src_pattern = os.path.join(
                self.job_dir,
                os.path.join(*output_record['location']),
                output_record['basename']
            ).format(runfolder=self.sample_id, sample_id=self.user_sample_id)
            source_files = glob(src_pattern)

            # TODO: turn off renaming for multiple files, md5 checksum
            for f in source_files:
                dest = os.path.join(
                    output_loc,
                    output_record.get('new_name', os.path.basename(f))
                ).format(sample_id=self.user_sample_id)
                exit_status += util.transfer_output_file(f, dest)

        return exit_status


class ProcessSamples(Stage):
    lf_name = '.process_samples.done'
    sample_sheet = luigi.Parameter()

    def requires(self):
        for sample_project in self.sample_sheet.sample_projects:
            for sample_id in self.sample_sheet.sample_projects[sample_project].sample_ids:
                yield DataOutput(sample_project=sample_project, sample_id=sample_id)

    def execute(self):
        for req in self.requires():
            exit_status = yield req
            app_logger.info('Finished transfer of %s with exit status %s' % (req.sample_id, exit_status))

    def output(self):
        return self.input()


def pipeline(input_run_folder):
    reader.transform_sample_sheet(input_run_folder)
    sample_sheet = reader.SampleSheet(input_run_folder)

    luigi.run(
        cmdline_args=[
            '--input-run-folder', input_run_folder,
            '--Stage-input-run-folder', input_run_folder,
            '--Stage-job-dir', os.path.join(cfg['jobs_dir'], os.path.basename(input_run_folder)),
            '--DataOutput-sample-sheet', sample_sheet
        ],
        main_task_cls=DataOutput,
        local_scheduler=True
    )
