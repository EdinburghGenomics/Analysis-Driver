import shutil
from time import sleep
from os.path import join, exists
from egcg_core import executor, clarity, util
from analysis_driver import reader
from analysis_driver.util import bash_commands
from analysis_driver.config import default as cfg
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.exceptions import PipelineError, SequencingRunError
from analysis_driver.report_generation.report_crawlers import RunCrawler
from analysis_driver.transfer_data import output_run_data
from . import Stage

app_logger = log_cfg.get_logger('fastq_production')


class FqProductionStage(Stage):
    def _run(self):
        raise NotImplementedError

    @property
    def fastq_dir(self):
        return join(self.job_dir, 'fastq')

    @property
    def fastqs(self):
        return util.find_all_fastqs(self.fastq_dir)

    @property
    def fastq_pairs(self):
        return util.find_all_fastq_pairs(self.fastq_dir)


class Setup(FqProductionStage):
    def _run(self):
        run_info = reader.RunInfo(self.input_dir)
        reader.transform_sample_sheet(self.input_dir, remove_barcode=not run_info.mask.has_barcodes)
        sample_sheet = reader.SampleSheet(join(self.input_dir, 'SampleSheet_analysis_driver.csv'), has_barcode=run_info.mask.has_barcodes)
        if not sample_sheet.validate(run_info.mask):
            raise PipelineError('Validation failed. Check SampleSheet.csv and RunInfo.xml.')

        self.cache_data('sample_sheet', sample_sheet)
        self.cache_data('run_info', run_info)

        crawler = RunCrawler(self.dataset_name, sample_sheet)
        crawler.send_data()
        sleep(cfg.get('lims_delay', 900))

        run_status = clarity.get_run(self.dataset_name).udf.get('Run Status')
        # TODO: catch bcl2fastq error logs instead
        if run_status != 'RunCompleted':
            app_logger.error('Run status is \'%s\'. Stopping.', run_status)
            raise SequencingRunError(run_status)
        return 0


class Bcl2Fastq(FqProductionStage):
    previous_stages = (Setup,)

    def _run(self):
        run_info = self.get_cached_data('run_info')
        sample_sheet = self.get_cached_data('sample_sheet')

        mask = sample_sheet.generate_mask(run_info.mask)
        app_logger.info('bcl2fastq mask: ' + mask)  # e.g: mask = 'y150n,i6,y150n'

        return executor.execute(
            bash_commands.bcl2fastq(
                self.input_dir,
                self.fastq_dir,
                sample_sheet.filename,
                mask
            ),
            job_name='bcl2fastq',
            working_dir=self.job_dir,
            cpus=8,
            mem=32
        ).join()


class WellDups(FqProductionStage):
    previous_stages = (Bcl2Fastq,)

    def _run(self):
        # this could be executed at the same time as bcl2fastq, but we need the fastq directory to exist
        output_file = join(self.fastq_dir, self.dataset_name + '.wellduplicate')
        output_err = join(self.fastq_dir, self.dataset_name + '.wellduplicate.err')
        coord_file = cfg.query('well_duplicate', 'coord_file')

        cmd = cfg.query('tools', 'well_duplicate') + ' -f %s -r %s -s hiseq_x > %s 2> %s' % (
            coord_file,
            self.input_dir,
            output_file,
            output_err
        )
        return executor.execute(
            cmd,
            job_name='welldup',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


class QCStage(FqProductionStage):
    previous_stages = (Bcl2Fastq,)

    def _run(self):
        raise NotImplementedError


class SickleFilter(QCStage):
    def _run(self):
        # filter the adapter dimer from fastq with sickle
        return executor.execute(
            *[bash_commands.sickle_paired_end_in_place(fqs) for fqs in self.fastq_pairs],
            job_name='sickle_filter',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class Fastqc(QCStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fq) for fq in self.fastqs],
            job_name='fastqc',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class FqChk(QCStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.seqtk_fqchk(fq) for fq in self.fastqs],
            job_name='fqchk',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


class MD5Sum(QCStage):
    def _run(self):
        return executor.execute(
            *[bash_commands.md5sum(fq) for fq in self.fastqs],
            job_name='md5sum',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


class OutputQCData(FqProductionStage):
    previous_stages = (SickleFilter, Fastqc, FqChk, MD5Sum)  # , WellDups)

    def _run(self):
        # copy the samplesheet, runinfo and run_parameters.xml to the fastq dir
        for f in ('SampleSheet.csv', 'SampleSheet_analysis_driver.csv',
                  'runParameters.xml', 'RunInfo.xml', 'RTAConfiguration.xml'):
            shutil.copy2(join(self.input_dir, f), join(self.fastq_dir, f))
        if not exists(join(self.fastq_dir, 'InterOp')):
            shutil.copytree(join(self.input_dir, 'InterOp'), join(self.fastq_dir, 'InterOp'))

        # find conversion xml file and send the results to the Rest API
        conversion_xml = join(self.fastq_dir, 'Stats', 'ConversionStats.xml')
        if exists(conversion_xml):
            app_logger.info('Found ConversionStats. Sending data.')
            crawler = RunCrawler(self.dataset_name, self.get_cached_data('sample_sheet'), conversion_xml, self.fastq_dir)
            # TODO: review whether we need this
            json_file = join(self.fastq_dir, 'demultiplexing_results.json')
            crawler.write_json(json_file)
            crawler.send_data()
            return 0
        else:
            app_logger.error('File not found: %s', conversion_xml)
            return 1


class DataOutput(FqProductionStage):
    previous_stages = (OutputQCData,)

    def _run(self):
        return output_run_data(self.fastq_dir, self.dataset_name)
