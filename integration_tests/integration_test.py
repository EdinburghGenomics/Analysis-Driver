import os
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
from integration_tests.mocked_data import patch_pipeline

cfg.load_config_file(os.getenv('ANALYSISDRIVERCONFIG'), env_var='ANALYSISDRIVERENV')
integration_cfg = Configuration(os.getenv('INTEGRATIONCONFIG'))


def _now():
    return datetime.utcnow().strftime('%d/%m/%Y %H:%M:%S')


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
        rest_communication.post_entry(
            'projects',
            {'project_id': '10015AT', 'samples': ['10015AT0004']}
        )

        # clean up any previous tests
        self._try_rm_dir(os.path.join(cfg['run']['output_dir'], run_id))
        self._try_rm_dir(os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004'))
        self._try_rm_dir(os.path.join(cfg['jobs_dir'], run_id))
        self._try_rm_dir(os.path.join(cfg['jobs_dir'], '10015AT0004'))

        self._test_success = True

    def tearDown(self):
        assert self.container_id
        self._execute('docker', 'stop', self.container_id)
        self._execute('docker', 'rm', self.container_id)
        logging_default.handlers = set()
        logging_default.loggers = {}

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
        d = {s['stage_name']: s['exit_status'] for s in stages}
        self.expect_equal(d, {s: 0 for s in stage_names}, 'stages')

    @staticmethod
    def _try_rm_dir(path):
        if os.path.isdir(path):
            rmtree(path)

    @classmethod
    def _check_md5(cls, fp):
        if os.path.isfile(fp):
            if fp.endswith('.gz'):
                return cls._execute('zcat ' + fp + ' | md5sum', shell=True).split()[0]
            elif fp.endswith('.bam'):
                return cls._execute('samtools view -h ' + fp + ' | md5sum', shell=True).split()[0]
            elif os.path.isfile(fp + '.md5'):
                with open(fp + '.md5', 'r') as f:
                    return f.readline().split(' ')[0]
            else:
                return cls._execute('md5sum', fp).split()[0]

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

    @staticmethod
    def _execute(*cmd, shell=False):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell)
        out, err = p.communicate()
        if err:
            raise ValueError(err)
        return out.decode('utf-8').rstrip('\n')

    def test_demultiplexing(self):
        with patch_pipeline():
            exit_status = client.main(['--run'])
            self.assertEqual(exit_status, 0)

            self.expect_equal(
                sorted(rest_communication.get_document('projects')['samples']),
                ['10015AT000' + str(i) for i in (1, 2, 3, 4, 6, 7, 8, 9)],
                'project_samples'
            )
            self.expect_equal(
                len(rest_communication.get_document('samples', where={'sample_id': '10015AT0004'})['run_elements']),
                8
            )
            output_dir = os.path.join(cfg['run']['output_dir'], '150723_E00306_0025_BHCHK3CCXX')
            output_fastqs = util.find_all_fastqs(output_dir)
            self.expect_equal(len(output_fastqs), 126, '# fastqs')  # 14 undetermined + 112 samples
            self.expect_output_files(
                integration_cfg['demultiplexing']['files'],
                base_dir=output_dir
        )
            self.expect_qc_data(
                rest_communication.get_document(
                    'run_elements',
                    where={'run_element_id': '150723_E00306_0025_BHCHK3CCXX_1_GAGATTCC'}
                ),
                integration_cfg['demultiplexing']['qc']
            )
            self.expect_stage_data(integration_cfg['demultiplexing']['stages'])

        assert self._test_success

    def test_bcbio(self):
        with patch_pipeline():
            exit_status = client.main(['--sample'])
            self.assertEqual(exit_status, 0)

            # Rest data
            self.expect_qc_data(
                rest_communication.get_document('samples', where={'sample_id': '10015AT0004'}),
                integration_cfg['bcbio']['qc']
            )

            # md5s
            self.expect_output_files(
                integration_cfg['bcbio']['files'],
                base_dir=os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004')
            )

            self.expect_stage_data(integration_cfg['bcbio']['stages'])

        assert self._test_success

    def test_var_calling(self):
        with patch_pipeline(species='Canis lupus familiaris', analysis_type='Variant Calling'):
            exit_status = client.main(['--sample'])
            self.assertEqual(exit_status, 0)

            self.expect_qc_data(
                rest_communication.get_document('samples', where={'sample_id': '10015AT0004'}),
                integration_cfg['var_calling']['qc']
            )

            self.expect_output_files(
                integration_cfg['var_calling']['files'],
                base_dir=os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004')
            )
            self.expect_stage_data(integration_cfg['var_calling']['stages'])

        assert self._test_success

    def test_qc(self):
        with patch_pipeline(species='Canis lupus familiaris', analysis_type='Not Variant Calling'):
            exit_status = client.main(['--sample'])
            self.assertEqual(exit_status, 0)

            self.expect_qc_data(
                rest_communication.get_document('samples', where={'sample_id': '10015AT0004'}),
                integration_cfg['qc']['qc']
            )

            self.expect_output_files(
                integration_cfg['qc']['files'],
                base_dir=os.path.join(cfg['sample']['output_dir'], '10015AT', '10015AT0004')
            )
            self.expect_stage_data(integration_cfg['qc']['stages'])

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
