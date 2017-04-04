import os
import sys
import gzip
import json
import hashlib
import pytest
import subprocess
import argparse
from io import StringIO
from time import sleep
from shutil import rmtree
from datetime import datetime
from contextlib import redirect_stdout, contextmanager
from egcg_core.config import cfg, Configuration
from egcg_core.app_logging import logging_default
from egcg_core import rest_communication, notifications, util, executor
from unittest import TestCase
from unittest.mock import Mock, patch
from analysis_driver import client

cfg.load_config_file(os.getenv('ANALYSISDRIVERCONFIG'), env_var='ANALYSISDRIVERENV')
integration_cfg = Configuration(os.getenv('INTEGRATIONCONFIG'))
entry_point = sys.argv[0]


def _now():
    return datetime.utcnow().strftime('%d/%m/%Y %H:%M:%S')


def _fake_get_list_of_samples(sample_names):
    samples = []
    for n in sample_names:
        m = Mock(udf={'Yield for Quoted Coverage (Gb)': 0.9})
        m.name = n
        samples.append(m)
    return samples


def _fake_get_user_sample_id(sample_name, lenient=False):
    return 'uid_' + sample_name


def _fake_get_plate_id_and_well(sample_name):
    return [sample_name + '_plate', 1337]


def _fake_welldups(self):
    output_file = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate')
    output_err = os.path.join(self.output_directory, self.dataset.name + '.wellduplicate.err')
    coord_file = cfg.query('well_duplicate', 'coord_file')

    cmd = cfg.query('tools', 'well_duplicate') + ' -f %s -r %s -t 1101 -s hiseq_x > %s 2> %s' % (
        coord_file, self.run_directory, output_file, output_err
    )
    return executor.execute(cmd, job_name='welldup', working_dir=self.working_dir, cpus=1, mem=2, log_commands=False).join()


@contextmanager
def patch_pipeline(species='Homo sapiens', analysis_type='Variant Calling gatk'):
    patches = []

    def _patch(ppath, **kwargs):
        p = patch('analysis_driver.' + ppath, **kwargs)
        p.start()
        patches.append(p)

    def _fake_get_sample(sample_name):
        return Mock(name=sample_name, udf={'Coverage': 1337, 'Analysis Type': analysis_type})

    _patch('pipelines.clarity.get_species_from_sample', return_value=species)
    _patch('pipelines.clarity.get_sample', new=_fake_get_sample)
    _patch('report_generation.report_crawlers.clarity.get_species_from_sample', return_value=species)
    _patch('report_generation.report_crawlers.clarity.get_sample', new=_fake_get_sample)
    _patch('dataset.clarity.get_expected_yield_for_sample', return_value=0.9)
    _patch('report_generation.report_crawlers.clarity.get_expected_yield_for_sample', return_value=0.9)
    _patch('dataset_scanner.get_list_of_samples', new=_fake_get_list_of_samples)
    _patch('pipelines.demultiplexing.clarity.get_run', return_value=Mock(udf={'Run Status': 'RunCompleted'}))
    _patch('pipelines.common.clarity.find_project_name_from_sample', return_value='10015AT')
    _patch('pipelines.common.clarity.get_sample', return_value=Mock(udf={}))
    _patch('quality_control.genotype_validation.clarity.find_project_name_from_sample', return_value='10015AT')
    _patch('pipelines.bcbio_pipelines.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('pipelines.common.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('pipelines.qc_pipelines.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('pipelines.variant_calling.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('report_generation.report_crawlers.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('report_generation.report_crawlers.clarity.get_plate_id_and_well', new=_fake_get_plate_id_and_well)
    _patch('report_generation.report_crawlers.clarity.get_sample_gender')
    _patch('quality_control.genotype_validation.clarity.get_samples_arrived_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_samples_genotyped_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_samples_sequenced_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_sample_names_from_project', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_sample_genotype', return_value=set())
    _patch('quality_control.lane_duplicates.WellDuplicates._well_duplicates', new=_fake_welldups)
    _patch('pipelines.demultiplexing.time.sleep')

    yield

    for p in patches:
        p.stop()


class IntegrationTest(TestCase):
    container_id = None

    def setUp(self):
        assert self.container_id is None
        self.container_id = self._execute('docker', 'run', '-d', 'egcg_reporting_app')
        assert self.container_id
        container_info = json.loads(self._execute('docker', 'inspect', self.container_id))[0]
        # for now, assume the container is running on the main 'bridge' network
        container_ip = container_info['NetworkSettings']['Networks']['bridge']['IPAddress']
        container_port = list(container_info['Config']['ExposedPorts'])[0].rstrip('/tcp')
        container_url = 'http://' + container_ip + ':' + container_port + '/api/0.1'
        rest_communication.default._baseurl = container_url

        sleep(15)  # allow time for the container's database and API to start running

        run_id = '150723_E00306_0025_BHCHK3CCXX'
        barcode = 'GAGATTCC'

        run_elements = []
        for lane in range(1, 9):
            run_elements.append(
                {'run_id': run_id, 'project_id': '10015AT', 'sample_id': '10015AT0004',
                 'library_id': 'LP6002014-DTP_A04', 'run_element_id': '%s_%s_%s' % (run_id, lane, barcode),
                 'useable': 'yes', 'barcode': barcode, 'lane': lane, 'clean_q30_bases_r1': 57000000,
                 'clean_q30_bases_r2': 57000000, 'clean_reads': 1}
            )
        for e in run_elements:
            rest_communication.post_entry('run_elements', e)

        rest_communication.post_entry(
            'samples',
            {'library_id': 'LP6002014-DTP_A04', 'project_id': '10015AT', 'sample_id': '10015AT0004',
             'run_elements': [e['run_element_id'] for e in run_elements]}
        )

        # clean up any previous tests
        self._try_rm_dir(os.path.join(cfg['run']['output_dir'], run_id))
        self._try_rm_dir(os.path.join(cfg['sample']['output_dir'], '10015AT0004'))
        self._try_rm_dir(os.path.join(cfg['jobs_dir'], run_id))
        self._try_rm_dir(os.path.join(cfg['jobs_dir'], '10015AT0004'))

        self._test_success = True

    def tearDown(self):
        assert self.container_id
        self._execute('docker', 'stop', self.container_id)
        self._execute('docker', 'rm', self.container_id)
        logging_default.handlers = set()
        logging_default.loggers = {}

    def expect_equal(self, obs, exp, name=None):
        if obs != exp:
            print(
                'Equality check {}failed:\nobs: {}\nexp: {}'.format(
                    str(name) + ' ' if name else '',  str(obs), str(exp)
                )
            )
            self._test_success = False

    @staticmethod
    def _try_rm_dir(path):
        if os.path.isdir(path):
            rmtree(path)

    @staticmethod
    def _read_md5_file(md5_file):
        if os.path.isfile(md5_file):
            with open(md5_file, 'r') as f:
                return f.readline().split(' ')[0]

    @staticmethod
    def _execute(*cmd):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        if err:
            raise ValueError(err)
        return out.decode('utf-8').rstrip('\n')

    def test_demultiplexing(self):
        with patch_pipeline():
            sys.argv = [entry_point, '--run']
            exit_status = client.main()
            self.expect_equal(exit_status, 0)

            rest_communication.patch_entries('run_elements', {'useable': 'yes'}, where={'sample_id': '10015AT0004'})
            self.expect_equal(
                sorted(rest_communication.get_document('projects')['samples']),
                ['10015AT000' + str(i) for i in (1, 2, 3, 4, 6, 7, 8, 9)],
                'project_samples'
            )
            self.expect_equal(
                len(rest_communication.get_document('samples', where={'sample_id': '10015AT0004'})['run_elements']),
                8,
                'run_element_count'
            )
            output_dir = os.path.join(cfg['run']['output_dir'], '150723_E00306_0025_BHCHK3CCXX')
            output_fastqs = util.find_files(output_dir, '*.fastq.gz') + util.find_files(output_dir, '10015AT', '*', '*.fastq.gz')
            self.expect_equal(len(output_fastqs), 126, 'fastqs')  # 14 undetermined + 112 sample

            exp_md5s = [integration_cfg['demultiplexing']['md5s'].get(os.path.basename(f)) for f in output_fastqs]
            obs_md5s = [hashlib.md5(gzip.open(f, 'r').read()).hexdigest() for f in output_fastqs]
            self.expect_equal(obs_md5s, exp_md5s, 'md5s')
        assert self._test_success

    def test_bcbio(self):
        with patch_pipeline():
            sys.argv = [entry_point, '--sample']
            exit_status = client.main()
            self.expect_equal(exit_status, 0)

            # Rest data
            e = rest_communication.get_document('samples', where={'sample_id': '10015AT0004'})
            obs = (e['bam_file_reads'], e['coverage']['genome_size'], e['duplicate_reads'], e['called_gender'])
            qc = integration_cfg['bcbio']['qc']
            exp = (qc['bam_file_reads'], qc['genome_size'], qc['duplicate_reads'], qc['called_gender'])
            self.expect_equal(obs, exp, 'qc')

            # md5s
            output_dir = os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004')
            files = ('samtools_stats.txt', 'uid_10015AT0004.vcf.stats', 'programs.txt',
                     'project-summary.yaml', 'uid_10015AT0004_R1_screen.txt', 'taxa_identified.json')
            obs = [self._read_md5_file(os.path.join(output_dir, f + '.md5')) for f in files]
            exp = [integration_cfg['bcbio']['md5s'].get(f) for f in files]
            self.expect_equal(obs, exp, 'md5s')
        assert self._test_success

    def test_var_calling(self):
        with patch_pipeline(species='Canis lupus familiaris', analysis_type='Variant Calling'):
            sys.argv = [entry_point, '--sample']
            exit_status = client.main()
            self.expect_equal(exit_status, 0)
            obs_sample = rest_communication.get_document('samples', where={'sample_id': '10015AT0004'})

            qcs = ('coverage', 'expected_yield', 'species_name', 'duplicate_reads', 'mapped_reads',
                   'properly_mapped_reads', 'bam_file_reads', 'median_coverage', 'species_contamination')
            obs = [obs_sample.get(k) for k in qcs]
            exp = [integration_cfg['var_calling']['qc'].get(k) for k in qcs]
            self.expect_equal(obs, exp, 'qc')

            output_dir = os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004')
            files = ('samtools_stats.txt', 'uid_10015AT0004.depth', 'uid_10015AT0004_R1_screen.txt',
                     'taxa_identified.json')
            obs = [self._read_md5_file(os.path.join(output_dir, f + '.md5')) for f in files]
            exp = [integration_cfg['var_calling']['md5s'].get(f) for f in files]
            self.expect_equal(obs, exp, 'md5s')
        assert self._test_success

    def test_qc(self):
        with patch_pipeline(species='Canis lupus familiaris', analysis_type='Not Variant Calling'):
            sys.argv = [entry_point, '--sample']
            exit_status = client.main()
            self.expect_equal(exit_status, 0)
            obs_sample = rest_communication.get_document('samples', where={'sample_id': '10015AT0004'})

            qcs = ('coverage', 'expected_yield', 'species_name', 'duplicate_reads', 'mapped_reads',
                   'properly_mapped_reads', 'bam_file_reads', 'median_coverage', 'species_contamination')
            obs = [obs_sample.get(k) for k in qcs]
            exp = [integration_cfg['var_calling']['qc'].get(k) for k in qcs]
            self.expect_equal(obs, exp, 'qc')

            output_dir = os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004')
            files = ('samtools_stats.txt', 'uid_10015AT0004.depth', 'uid_10015AT0004_R1_screen.txt',
                     'taxa_identified.json')
            obs = [self._read_md5_file(os.path.join(output_dir, f + '.md5')) for f in files]
            exp = [integration_cfg['var_calling']['md5s'].get(f) for f in files]
            self.expect_equal(obs, exp, 'md5s')
        assert self._test_success


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--noemail', dest='email', action='store_false')
    args = a.parse_args()

    start_time = _now()
    s = StringIO()
    with redirect_stdout(s):
        pytest.main([__file__])
    end_time = _now()

    test_output = util.str_join(
        'Pipeline test finished. ',
        'Start time: %s, finish time: %s. '
        'Pytest output:\n' % (start_time, end_time),
        s.getvalue()
    )
    print(test_output)
    if args.email:
        e = notifications.EmailNotification(
            'Analysis Driver integration test',
            **integration_cfg['notification']
        )
        e.notify(test_output)


if __name__ == '__main__':
    main()
