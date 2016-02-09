__author__ = 'mwham'
import luigi
import os.path
import shutil
import yaml
from glob import glob
from analysis_driver import executor, writer, reader, util, clarity
from analysis_driver.report_generation import report_crawlers
from analysis_driver.transfer_data import prepare_run_data, prepare_sample_data
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('luigi')


class Stage(luigi.Task):
    dataset = luigi.Parameter()

    @property
    def job_dir(self):
        return os.path.join(cfg['jobs_dir'], self.dataset.name)

    @property
    def input_dir(self):
        return os.path.join(cfg['input_dir'], self.dataset.name)

    def output(self):
        return luigi.LocalTarget(self._stage_lock_file())

    def run(self):
        exit_status = self._run()
        if exit_status == 0:
            self._touch(self._stage_lock_file())
        else:
            raise AnalysisDriverError('Exit status was %s. Stopping' % exit_status)

    def _run(self):
        raise NotImplementedError

    @staticmethod
    def _touch(f):
        open(f, 'w').close()

    def _stage_lock_file(self):
        return os.path.join(self.job_dir, '.' + self.__class__.__name__.lower() + '.done')


class TransferRun(Stage):
    def _run(self):
        dataset_dir = prepare_run_data(dataset=self.dataset)
        if os.path.isdir(dataset_dir):
            return 0
        else:
            return 1  # TODO: return an exit status from prepare_run_data instead


class Bcl2Fastq(Stage):
    def _run(self):
        fastq_dir = os.path.join(self.job_dir, 'fastq')

        sample_sheet = reader.SampleSheet(os.path.join(self.input_dir, 'SampleSheet_analysis_driver.csv'))
        run_info = reader.RunInfo(self.input_dir)
        mask = sample_sheet.generate_mask(run_info.mask)
        bcl2fastq = writer.bash_commands.bcl2fastq(
            self.input_dir,
            fastq_dir,
            sample_sheet.filename,
            mask
        )

        return executor.execute(
            [bcl2fastq],
            job_name='bcl2fastq',
            run_id=self.dataset.name,
            cpus=8,
            mem=32
        ).join()

    @staticmethod
    def requires():
        return TransferRun()


class DemultiplexingFastqc(Stage):
    def _run(self):
        fastq_dir = os.path.join(self.job_dir, 'fastq')
        return executor.execute(
            [writer.bash_commands.fastqc(fq) for fq in util.fastq_handler.find_all_fastqs(fastq_dir)],
            job_name='fastqc',
            run_id=self.dataset.name,
            cpus=1,
            mem=2
        ).join()

    @staticmethod
    def requires():
        return Bcl2Fastq()


class DemultiplexingMD5Sum(Stage):
    def _run(self):
        fastq_dir = os.path.join(self.job_dir, 'fastq')

        return executor.execute(
            [writer.bash_commands.md5sum(fq) for fq in util.fastq_handler.find_all_fastqs(fastq_dir)],
            job_name='md5sum',
            run_id=self.dataset.name,
            cpus=1,
            mem=2,
            log_command=False
        ).join()

    @staticmethod
    def requires():
        return Bcl2Fastq()


class DemultiplexingOutput(Stage):
    def _run(self):
        exit_status = 0
        fastq_dir = os.path.join(self.job_dir, 'fastq')

        transfer_files = (
            'SampleSheet.csv', 'SampleSheet_analysis_driver.csv', 'runParameters.xml', 'RunInfo.xml',
            'RTAConfiguration.xml'
        )
        for f in transfer_files:
            shutil.copy2(os.path.join(self.input_dir, f), os.path.join(fastq_dir, f))
        if not os.path.exists(os.path.join(fastq_dir, 'InterOp')):
            shutil.copytree(os.path.join(self.input_dir, 'InterOp'), os.path.join(fastq_dir, 'InterOp'))

        conversion_xml = os.path.join(fastq_dir, 'Stats', 'ConversionStats.xml')
        if os.path.exists(conversion_xml):
            app_logger.info('Found ConversionStats. Sending data.')
            crawler = report_crawlers.RunCrawler(
                self.dataset.name,
                reader.SampleSheet(os.path.join(self.input_dir, 'SampleSheet_analysis_driver.csv')),
                conversion_xml
            )
            # TODO: review whether we need this
            json_file = os.path.join(fastq_dir, 'demultiplexing_results.json')
            crawler.write_json(json_file)
            crawler.send_data()
        else:
            app_logger.error('File not found: %s' % conversion_xml)
            exit_status += 1

        return exit_status

    @staticmethod
    def requires():
        return DemultiplexingFastqc(), DemultiplexingMD5Sum()


class BCBioPrepareSamples(Stage):
    def _run(self):
        fastqs = prepare_sample_data(self.dataset.name)
        user_sample_id = clarity.get_user_sample_name(self.dataset.name, lenient=True)
        cmd = util.bcbio_prepare_samples_cmd(
            self.job_dir,
            self.dataset.name,
            fastqs,
            user_sample_id
        )
        return executor.execute([cmd], job_name='bcbio_prepare_samples', run_id=self.dataset.name).join()


class MergedFastqc(Stage):
    def _run(self):
        user_sample_id = clarity.get_user_sample_name(self.dataset.name, lenient=True)
        fastq_pair = glob(os.path.join(self.job_dir, 'merged', user_sample_id + '_R?.fastq.gz'))
        return executor.execute(
            [writer.bash_commands.fastqc(fastq_file) for fastq_file in fastq_pair],
            job_name='fastqc2',
            run_id=self.dataset.name,
            cpus=1,
            mem=2
        )

    @staticmethod
    def requires():
        return BCBioPrepareSamples()



class BCBio(Stage):
    sample_project = luigi.Parameter()
    sample_id = luigi.Parameter()
    fastqs = luigi.Parameter()

    @staticmethod
    def requires():
        return BCBioPrepareSamples()

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

    def modify_run_yaml(self):
        bcbio_dir = os.path.join(self.job_dir, 'samples_' + self.dataset.name + '-merged')
        run_yaml = os.path.join(bcbio_dir, 'config', 'samples_' + self.dataset.name + '-merged.yaml')

        user_sample_id = clarity.get_user_sample_name(self.dataset.name)
        if user_sample_id:
            app_logger.debug('Found user sample: ' + user_sample_id)

            with open(run_yaml, 'r') as i:
                run_config = yaml.load(i)
            run_config['fc_name'] = user_sample_id
            with open(run_yaml, 'w') as o:
                o.write(yaml.safe_dump(run_config, default_flow_style=False))
        return 0

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
            prelim_cmds=writer.bash_commands.export_env_vars(),
            job_name='bcbio',
            run_id=self.run_id,
            cpus=10,
            mem=64
        ).join()

    def _run(self):
        for stage in (self.prepare_template, self.modify_run_yaml, self.run_bcbio):
            ex = stage()
            if ex != 0:
                return ex
        return 0


class BCBioOutput(Stage):

    @staticmethod
    def requires():
        return BCBio(), MergedFastqc()

    def _run(self):

        output_dir = os.path.join(
            cfg['output_dir'],
            clarity.find_project_from_sample(self.dataset.name),
            self.dataset.name
        )
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        bcbio_source_dir = os.path.join(
            self.job_dir,
            'samples_' + self.dataset.name + '-merged',
            'final',
            self.dataset.name
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
        return util.transfer_output_files(self.dataset.name, output_dir, source_path_mapping)


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
