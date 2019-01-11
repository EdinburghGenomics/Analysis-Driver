import os
import subprocess
from egcg_core import rest_communication, util
from egcg_core.config import cfg
from egcg_core.app_logging import logging_default
from egcg_core.integration_testing import ReportingAppIntegrationTest
from unittest.mock import Mock, patch
from analysis_driver import client
from integration_tests import mocked_data
from integration_tests.mocked_data import MockedRunProcess, mocked_flowcell_pooling


class IntegrationTest(ReportingAppIntegrationTest):
    patches = (
        patch('analysis_driver.client.load_config'),
        patch('egcg_core.clarity.find_project_name_from_sample', return_value='10015AT'),
        patch('egcg_core.clarity.get_plate_id_and_well', new=mocked_data.fake_get_plate_id_and_well),
        patch('egcg_core.clarity.get_project', return_value=mocked_data.mocked_clarity_project),
        patch('egcg_core.clarity.get_run', return_value=mocked_data.mocked_clarity_run),
        patch('egcg_core.clarity.get_sample_gender'),
        patch('egcg_core.clarity.get_sample_genotype', return_value=set()),
        patch('egcg_core.clarity.get_sample_names_from_project', return_value=set()),
        patch('egcg_core.clarity.get_samples_arrived_with', return_value=set()),
        patch('egcg_core.clarity.get_samples_genotyped_with', return_value=set()),
        patch('egcg_core.clarity.get_samples_sequenced_with', return_value=set()),
        patch('egcg_core.clarity.get_user_sample_name', new=mocked_data.fake_get_user_sample_id),
        patch('analysis_driver.dataset.LimsNotification'),
        patch(
            'analysis_driver.quality_control.well_duplicates.WellDuplicates._welldups_cmd',
            new=mocked_data.fake_welldups_cmd
        ),
        patch(
            'analysis_driver.pipelines.demultiplexing.BadTileCycleDetector',
            return_value=Mock(
                detect_bad_cycles=Mock(return_value={5: [309, 310]}),
                detect_bad_tiles=Mock(return_value={})
            )
        )
    )

    def __init__(self, *args):
        super().__init__(*args)
        self.run_id = self.cfg['input_data']['run_id']
        self.barcode = self.cfg['input_data']['barcode']
        self.project_id = self.cfg['input_data']['project_id']
        self.sample_id = self.cfg['input_data']['sample_id']
        self.library_id = self.cfg['input_data']['library_id']
        self.samples_for_project = self.cfg['input_data']['samples_for_project']
        self.run_dir = os.path.dirname(os.getcwd())  # we're inside the checked out project, not the top level

    def setUp(self):
        super().setUp()
        cfg.load_config_file(os.getenv('ANALYSISDRIVERCONFIG'), env_var='ANALYSISDRIVERENV')
        run_elements = []
        for lane in range(1, 9):
            run_elements.append(
                {'run_id': self.run_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
                 'library_id': self.library_id, 'run_element_id': '%s_%s_%s' % (self.run_id, lane, self.barcode),
                 'useable': 'yes', 'barcode': self.barcode, 'lane': lane, 'bases_r1': 60000000, 'bases_r2': 60000000,
                 'clean_bases_r1': 58000000, 'clean_bases_r2': 58000000, 'q30_bases_r1': 57000000,
                 'q30_bases_r2': 57000000, 'clean_q30_bases_r1': 56000000, 'clean_q30_bases_r2': 56000000,
                 'clean_reads': 1}
            )
        for e in run_elements:
            rest_communication.post_entry('run_elements', e)

        rest_communication.post_entry(
            'samples',
            {'library_id': self.library_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
             'run_elements': [e['run_element_id'] for e in run_elements], 'required_yield': 900000000}
        )

        rest_communication.post_entry('projects', {'project_id': self.project_id, 'samples': [self.sample_id]})

        self.dynamic_patches = []
        self._test_success = True

    def tearDown(self):
        super().tearDown()
        self._reset_logging()
        for p in self.dynamic_patches:
            p.stop()

    def setup_test(self, test_type, test_name, integration_section, species='Homo sapiens', analysis_type='Variant Calling gatk'):
        cfg.content['jobs_dir'] = os.path.join(os.path.dirname(os.getcwd()), 'jobs', test_name)
        cfg.content[test_type]['output_dir'] = os.path.join(os.path.dirname(os.getcwd()), 'outputs', test_name)
        if 'input_dir' in self.cfg[integration_section]:
            cfg.content[test_type]['input_dir'] = self.cfg[integration_section]['input_dir']

        os.makedirs(cfg['jobs_dir'])
        os.makedirs(cfg[test_type]['output_dir'])

        def _fake_get_sample(sample_name):
            return Mock(
                name=sample_name,
                udf={
                    'Coverage': 1337,
                    'Analysis Type': analysis_type,
                    'Yield for Quoted Coverage (Gb)': 15,
                    'Required Yield (Gb)': 30,
                    'Coverage (X)': 15,
                }
            )

        for p in (patch('egcg_core.clarity.get_sample', new=_fake_get_sample),
                  patch('egcg_core.clarity.get_species_from_sample', return_value=species)):
            self.dynamic_patches.append(p)
            p.start()

    @staticmethod
    def _reset_logging():
        for logger in logging_default.loggers.values():
            for handler in logger.handlers:
                logger.removeHandler(handler)

        logging_default.handlers = set()
        logging_default.loggers = {}

    def expect_equal(self, obs, exp, name):
        try:
            self.assertEqual(name, obs, exp)
        except AssertionError:
            self._test_success = False

    def expect_output_files(self, exp, base_dir):
        for k, v in exp.items():
            if base_dir:
                k = os.path.join(base_dir, k)

            if v is None:
                self.expect_equal(os.path.isfile(k), True, k)
            else:
                self.expect_equal(self._check_md5(k), v, k)

    def expect_qc_data(self, obs, exp):
        for e in sorted(exp):
            self.expect_equal(util.query_dict(obs, e), exp[e], e)

    def expect_stage_data(self, stage_names, **query_kw):
        """
        Take a list of stage name as string assuming exit status is 0 or as tuple with (stage_name, exist status).
        Compare to the list of stage name and exist status retrieve from the analysis_driver_stages endpoint using the
        keyword arguments provided.
        """
        stages = rest_communication.get_documents('analysis_driver_stages', **query_kw)
        obs = {s['stage_name']: s.get('exit_status') for s in stages}
        exp = dict([s if type(s) == tuple else (s, 0) for s in stage_names])
        self.expect_equal(obs, exp, 'stages')

    def _check_md5(self, fp):
        if not os.path.isfile(fp):
            return None

        elif fp.endswith('.bam'):
            cmd = "{samtools} idxstats {fp} | awk '{awk_exp}'".format(
                samtools=self.cfg['samtools'],
                fp=fp,
                awk_exp='{seq_len+=$2; mapped+=$3; unmapped+=$4}END {print seq_len,mapped,unmapped}'
            )
        elif fp.endswith('.gz'):
            cmd = 'zcat ' + fp
        else:
            cmd = 'cat ' + fp

        cmd += " | sed 's/{cwd}//g' | md5sum".format(cwd=self.run_dir.replace('/', '\/'))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        out, err = p.communicate()

        if err:
            raise ValueError(err)

        return out.decode('utf-8').rstrip('\n').split()[0]

    def test_demultiplexing(self):
        self.setup_test('run', 'test_demultiplexing', 'demultiplexing')

        exit_status = client.main(['--run'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_equal(
            sorted(rest_communication.get_document('projects')['samples']),
            ['10015AT000' + str(i) for i in (1, 2, 3, 4, 6, 7, 8, 9)],
            'project samples'
        )
        self.expect_equal(
            len(rest_communication.get_document('samples', where={'sample_id': self.sample_id})['run_elements']),
            8,
            '# run elements'
        )
        output_dir = os.path.join(cfg['run']['output_dir'], self.run_id)
        output_fastqs = util.find_all_fastqs(output_dir)
        self.expect_equal(len(output_fastqs), self.cfg['demultiplexing']['nb_fastq'], '# fastqs')
        self.expect_output_files(
            self.cfg['demultiplexing']['files'],
            base_dir=output_dir
        )
        self.expect_qc_data(
            rest_communication.get_document(
                'run_elements',
                where={'run_element_id': self.run_id + '_1_' + self.barcode}
            ),
            self.cfg['demultiplexing']['run_element_qc']
        )
        if 'lane_qc' in self.cfg['demultiplexing']:
            self.expect_qc_data(
                rest_communication.get_document(
                    'lanes',
                    where={'lane_id': self.run_id + '_1'}
                ),
                self.cfg['demultiplexing']['lane_qc']
            )
        self.expect_stage_data(['setup', 'wellduplicates', 'bcl2fastq', 'phixdetection', 'fastqfilter', 'seqtkfqchk',
                                'md5sum', 'fastqc', 'integritycheck', 'qcoutput1', 'dataoutput', 'cleanup',
                                'samtoolsdepthmulti', 'picardinsertsizemulti', 'qcoutput2', 'runreview',
                                'picardmarkduplicatemulti', 'samtoolsstatsmulti', 'bwaalignmulti', 'picardgcbias'])

        proc = rest_communication.get_document('analysis_driver_procs')
        self.expect_equal(
            rest_communication.get_document('runs', where={'run_id': self.run_id}).get('analysis_driver_procs'),
            [proc['proc_id']],
            'run proc registered'
        )
        self.expect_equal(
            proc['pipeline_used'],
            {'name': 'demultiplexing', 'toolset_version': 0, 'toolset_type': 'run_processing'},
            'pipeline used'
        )
        assert self._test_success

    def test_demultiplexing_aborted(self):
        self.setup_test('sample', 'test_demultiplexing_aborted', 'demultiplexing')

        mocked_clarity_run = MockedRunProcess(udf={'Run Status': 'RunAborted'}, container=mocked_flowcell_pooling)
        with patch('egcg_core.clarity.get_run', return_value=mocked_clarity_run):
            exit_status = client.main(['--run'])

        self.assertEqual('exit status', exit_status, 0)
        ad_proc = rest_communication.get_document('analysis_driver_procs', where={'dataset_name': self.run_id})

        self.expect_equal(
            ad_proc['status'], 'aborted', 'pipeline status'
        )
        self.expect_stage_data([('setup', 9)], where={'analysis_driver_proc': ad_proc['proc_id']})

    def test_bcbio(self):
        self.setup_test('sample', 'test_bcbio', 'bcbio')
        exit_status = client.main(['--sample'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_qc_data(
            rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
            self.cfg['bcbio']['qc']
        )

        self.expect_output_files(
            self.cfg['bcbio']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        self.expect_stage_data(['mergefastqs', 'fastqc', 'genotypevalidation', 'bcbio', 'fastqscreen',
                                'fixunmapped', 'blast', 'gendervalidation', 'vcfstats', 'samtoolsdepth',
                                'verifybamid', 'sampledataoutput', 'md5sum', 'cleanup', 'samplereview'])

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['pipeline_used'],
            {'toolset_type': 'human_sample_processing', 'name': 'bcbio', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 9)],
            'data source'
        )

        assert self._test_success

    def test_var_calling(self):
        self.setup_test('sample', 'test_var_calling', 'var_calling', 'Canis lupus familiaris', 'Variant Calling')
        exit_status = client.main(['--sample'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_qc_data(
            rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
            self.cfg['var_calling']['qc']
        )

        self.expect_output_files(
            self.cfg['var_calling']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        self.expect_stage_data(['mergefastqs', 'samplereview', 'fastqscreen', 'printreads', 'blast', 'baserecal',
                                'samtoolsdepth', 'vcfstats', 'selectvariants', 'fastqc', 'variantfiltration',
                                'cleanup', 'realign', 'realigntarget', 'samtoolsstats', 'haplotypecaller',
                                'md5sum', 'bwamem', 'sampledataoutput', 'genotypegvcfs'])

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['pipeline_used'],
            {'toolset_type': 'non_human_sample_processing', 'name': 'variant_calling', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 9)],
            'data source'
        )

        assert self._test_success

    def _run_qc_test(self):
        exit_status = client.main(['--sample'])

        self.assertEqual('exit status', exit_status, 0)
        self.expect_qc_data(
            rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
            self.cfg['qc']['qc']
        )
        self.expect_stage_data(['vcfstats', 'selectvariants', 'fastqc', 'variantfiltration', 'samtoolsdepth',
                                'mergefastqs', 'samtoolsstats', 'samplereview', 'haplotypecaller', 'cleanup',
                                'fastqscreen', 'md5sum', 'bwamem', 'sampledataoutput', 'genotypegvcfs', 'blast'])

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['pipeline_used'],
            {'toolset_type': 'non_human_sample_processing', 'name': 'qc', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(
                rest_communication.get_document('analysis_driver_procs')['data_source'],
                ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 9)],
                'data source'
            )

    def test_qc(self):
        self.setup_test('sample', 'test_qc', 'qc', 'Canis lupus familiaris', 'Not Variant Calling')

        self._run_qc_test()
        self.expect_output_files(
            self.cfg['qc']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        assert self._test_success

    def test_resume(self):
        self.setup_test('sample', 'test_resume', 'qc', 'Canis lupus familiaris', 'Not Variant Calling')
        with patch('analysis_driver.pipelines.common.BWAMem._run', return_value=1):
            exit_status = client.main(['--sample'])
            self.assertEqual('exit status', exit_status, 9)

        procs = rest_communication.get_documents('analysis_driver_procs')
        self.expect_equal(len(procs), 1, '# procs')
        self.expect_equal(procs[0]['status'], 'failed', 'process failed')
        self.expect_equal(
            sorted(s['stage_name'] for s in rest_communication.get_documents('analysis_driver_stages')),
            ['bwamem', 'fastqc', 'mergefastqs'],
            'stages'
        )
        self.expect_equal(
            rest_communication.get_document('analysis_driver_stages', where={'stage_name': 'bwamem'})['exit_status'],
            1,
            'bwamem exit status'
        )

        self._reset_logging()
        client.main(['--sample', '--resume', '10015AT0004'])
        self.expect_equal(rest_communication.get_document('analysis_driver_procs')['status'], 'resume', 'resumed')
        self._reset_logging()

        self._run_qc_test()

        exp_files = self.cfg['qc']['files'].copy()
        exp_files.update(self.cfg['resume']['files'])
        self.expect_output_files(
            exp_files,
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        procs = rest_communication.get_documents('analysis_driver_procs')
        self.expect_equal(len(procs), 1, 'used existing proc')
        self.expect_equal(procs[0]['status'], 'finished', 'proc status finished')

        assert self._test_success

    def test_project(self):
        # samples for the project process tests
        for sample in self.samples_for_project:
            rest_communication.post_entry(
                'samples', {'project_id': self.project_id, 'sample_id': sample, 'required_yield': 120000000000,
                            'user_sample_id': 'uid_' + sample}
            )
            rest_communication.post_entry(
                'analysis_driver_procs',
                {'proc_id': sample + 'proc_id', 'dataset_name': sample, 'dataset_type': 'sample',
                 'status': 'finished', 'pipeline_used': {'name': 'bcbio'}}
            )
        rest_communication.patch_entry(
            'projects',
            {'samples': self.samples_for_project},
            'project_id',
            self.project_id,
            update_lists=['samples']
        )

        self.setup_test('project', 'test_project', 'project')
        exit_status = client.main(['--project'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_output_files(
            self.cfg['project']['files'],
            base_dir=os.path.join(cfg['project']['output_dir'], self.project_id)
        )

        self.expect_stage_data(['genotypegvcfs', 'relatedness', 'peddy', 'parserelatedness', 'md5sum', 'output',
                                'cleanup'])
        ad_procs = rest_communication.get_document('analysis_driver_procs', where={'dataset_name': self.project_id})
        self.expect_equal(
            ad_procs['pipeline_used'],
            {'toolset_type': 'project_processing', 'name': 'project', 'toolset_version': 0},
            'pipeline used'
        )

        assert self._test_success
