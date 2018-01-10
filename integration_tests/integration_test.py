import os
import sys
import json
import pytest
import subprocess
import argparse
from io import StringIO
from time import sleep
from shutil import rmtree
from datetime import datetime
from contextlib import redirect_stdout
from egcg_core import rest_communication, notifications, util
from egcg_core.config import cfg, Configuration
from egcg_core.app_logging import logging_default
from unittest import TestCase
from analysis_driver import client
from analysis_driver.config import load_config
from integration_tests.mocked_data import patch_pipeline

cfg.load_config_file(os.getenv('ANALYSISDRIVERCONFIG'), env_var='ANALYSISDRIVERENV')
integration_cfg = Configuration(os.getenv('INTEGRATIONCONFIG'))


def now():
    return datetime.utcnow().strftime('%Y-%m-%d_%H:%M:%S')


def execute(*cmd, shell=False):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell)
    out, err = p.communicate()
    if err:
        raise ValueError(err)
    return out.decode('utf-8').rstrip('\n')


class IntegrationTest(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.container_id = None

        cls.run_id = integration_cfg['input_data']['run_id']
        cls.barcode = integration_cfg['input_data']['barcode']
        cls.project_id = integration_cfg['input_data']['project_id']
        cls.sample_id = integration_cfg['input_data']['sample_id']
        cls.library_id = integration_cfg['input_data']['library_id']

        load_config()
        cls.original_job_dir = cfg['jobs_dir']
        cls.original_run_output = cfg['run']['output_dir']
        cls.original_sample_output = cfg['sample']['output_dir']

        # clean up any previous tests
        for top_level in (cls.original_job_dir, cls.original_run_output, cls.original_sample_output):
            for d in os.listdir(top_level):
                rmtree(os.path.join(top_level, d))

    def setUp(self):
        assert self.container_id is None
        self.container_id = execute(
            'docker', 'run', '-d', integration_cfg['reporting_app']['image_name'],
            integration_cfg.query('reporting_app', 'branch', ret_default='master')
        )
        assert self.container_id
        container_info = json.loads(execute('docker', 'inspect', self.container_id))[0]
        # for now, assume the container is running on the main 'bridge' network
        container_ip = container_info['NetworkSettings']['Networks']['bridge']['IPAddress']
        container_port = list(container_info['Config']['ExposedPorts'])[0].rstrip('/tcp')
        container_url = 'http://' + container_ip + ':' + container_port + '/api/0.1'
        rest_communication.default._baseurl = container_url

        sleep(30)  # allow time for the container's database and API to start running

        run_elements = []
        for lane in range(1, 9):
            run_elements.append(
                {'run_id': self.run_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
                 'library_id': self.library_id, 'run_element_id': '%s_%s_%s' % (self.run_id, lane, self.barcode),
                 'useable': 'yes', 'barcode': self.barcode, 'lane': lane, 'clean_q30_bases_r1': 57000000,
                 'clean_q30_bases_r2': 57000000, 'clean_reads': 1, 'clean_bases_r1': 60000000, 'clean_bases_r2': 60000000}
            )
        for e in run_elements:
            rest_communication.post_entry('run_elements', e)

        rest_communication.post_entry(
            'samples',
            {'library_id': self.library_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
             'run_elements': [e['run_element_id'] for e in run_elements], 'required_yield': 100000000}
        )
        rest_communication.post_entry('projects', {'project_id': self.project_id, 'samples': [self.sample_id]})

        self._test_success = True

    def tearDown(self):
        assert self.container_id
        execute('docker', 'stop', self.container_id)
        execute('docker', 'rm', self.container_id)

        for logger in logging_default.loggers.values():
            for handler in logger.handlers:
                logger.removeHandler(handler)

        logging_default.handlers = set()
        logging_default.loggers = {}

        cfg.content['jobs_dir'] = self.original_job_dir
        cfg.content['run']['output_dir'] = self.original_run_output
        cfg.content['sample']['output_dir'] = self.original_sample_output

    @staticmethod
    def setup_test(test_type, test_name):
        cfg.content['jobs_dir'] = os.path.join(cfg['jobs_dir'], test_name)
        cfg.content[test_type]['output_dir'] = os.path.join(cfg[test_type]['output_dir'], test_name)
        os.mkdir(cfg['jobs_dir'])
        os.mkdir(cfg[test_type]['output_dir'])

    def expect_equal(self, obs, exp, name=''):
        if name:
            name += ' '

        if obs != exp:
            print('Check %sfailed:\nobs: %s\nexp: %s' % (name, obs, exp))
            self._test_success = False

    def expect_output_files(self, exp, base_dir=None):
        for k, v in exp.items():
            if base_dir:
                k = os.path.join(base_dir, k)

            if v is None:
                self.expect_equal(os.path.isfile(k), True, k)
            else:
                self.expect_equal(self._check_md5(k), v, k)

    def expect_qc_data(self, obs, exp):
        self.expect_equal(
            [self._query_dict(obs, e) for e in sorted(exp)],
            [exp[e] for e in sorted(exp)]
        )

    def expect_stage_data(self, stage_names):
        stages = rest_communication.get_documents('analysis_driver_stages')
        obs = {s['stage_name']: s['exit_status'] for s in stages}
        self.expect_equal(obs, {s: 0 for s in stage_names}, 'stages')

    @staticmethod
    def _check_md5(fp):
        if os.path.isfile(fp):
            if fp.endswith('.gz'):
                return execute('zcat ' + fp + ' | md5sum', shell=True).split()[0]
            elif fp.endswith('.bam'):
                cmd = "{samtools} idxstats {fp} | awk '{awk_exp}' | md5sum".format(
                    samtools=integration_cfg['samtools'],
                    fp=fp,
                    awk_exp='{seq_len+=$2; mapped+=$3; unmapped+=$4}END {print seq_len,mapped,unmapped}'
                )
                return execute(cmd, shell=True).split()[0]
            elif os.path.isfile(fp + '.md5'):
                with open(fp + '.md5', 'r') as f:
                    return f.readline().split(' ')[0]
            else:
                return execute('md5sum', fp).split()[0]

    @staticmethod
    def _query_dict(input_dict, path):
        i = input_dict.copy()
        v = None
        for p in path.split('.'):
            v = i.get(p)
            if v is None:
                return None
            elif type(v) is list:
                i = {str(v.index(x)): x for x in v}
            elif type(v) is dict:
                i = v
        return v

    def test_demultiplexing(self):
        self.setup_test('run', 'test_demultiplexing')
        with patch_pipeline():
            exit_status = client.main(['--run'])
            self.assertEqual(exit_status, 0)

            self.expect_equal(
                sorted(rest_communication.get_document('projects')['samples']),
                ['10015AT000' + str(i) for i in (1, 2, 3, 4, 6, 7, 8, 9)],
                'project_samples'
            )
            self.expect_equal(
                len(rest_communication.get_document('samples', where={'sample_id': self.sample_id})['run_elements']),
                8
            )
            output_dir = os.path.join(cfg['run']['output_dir'], self.run_id)
            output_fastqs = util.find_all_fastqs(output_dir)
            self.expect_equal(len(output_fastqs), 126, '# fastqs')  # 14 undetermined + 112 samples
            self.expect_output_files(
                integration_cfg['demultiplexing']['files'],
                base_dir=output_dir
            )
            self.expect_qc_data(
                rest_communication.get_document(
                    'run_elements',
                    where={'run_element_id': self.run_id + '_1_' + self.barcode}
                ),
                integration_cfg['demultiplexing']['qc']
            )
            self.expect_stage_data(integration_cfg['demultiplexing']['stages'])

        assert self._test_success

    def test_bcbio(self):
        self.setup_test('sample', 'test_bcbio')
        with patch_pipeline():
            exit_status = client.main(['--sample'])
            self.assertEqual(exit_status, 0)

            self.expect_qc_data(
                rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
                integration_cfg['bcbio']['qc']
            )

            self.expect_output_files(
                integration_cfg['bcbio']['files'],
                base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
            )
            self.expect_stage_data(integration_cfg['bcbio']['stages'])

        assert self._test_success

    def test_var_calling(self):
        self.setup_test('sample', 'test_var_calling')
        with patch_pipeline(species='Canis lupus familiaris', analysis_type='Variant Calling'):
            exit_status = client.main(['--sample'])
            self.assertEqual(exit_status, 0)

            self.expect_qc_data(
                rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
                integration_cfg['var_calling']['qc']
            )

            self.expect_output_files(
                integration_cfg['var_calling']['files'],
                base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
            )
            self.expect_stage_data(integration_cfg['var_calling']['stages'])

        assert self._test_success

    def test_qc(self):
        self.setup_test('sample', 'test_qc')
        with patch_pipeline(species='Canis lupus familiaris', analysis_type='Not Variant Calling'):
            exit_status = client.main(['--sample'])
            self.assertEqual(exit_status, 0)

            self.expect_qc_data(
                rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
                integration_cfg['qc']['qc']
            )

            self.expect_output_files(
                integration_cfg['qc']['files'],
                base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
            )
            self.expect_stage_data(integration_cfg['qc']['stages'])

        assert self._test_success


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--stdout', action='store_true')
    a.add_argument('--email', action='store_true')
    a.add_argument('--log_repo')
    args = a.parse_args()

    start_time = now()
    s = StringIO()
    with redirect_stdout(s):
        exit_status = pytest.main([__file__])
    end_time = now()

    test_output = util.str_join(
        'Pipeline end-to-end test finished',
        'Run on commit %s' % execute('git', 'log', "--format=%h on%d, made on %aD", '-1'),
        'Start time: %s, finish time: %s' % (start_time, end_time),
        'Pytest output:',
        s.getvalue(),
        separator='\n'
    )

    if args.log_repo:
        with open(os.path.join(args.log_repo, start_time + '.log'), 'w') as f:
            f.write(test_output)

    if args.stdout:
        print(test_output)

    if args.email:
        notifications.send_email(
            test_output, subject='Analysis Driver integration test', **integration_cfg['notification']
        )

    return exit_status


if __name__ == '__main__':
    sys.exit(main())
