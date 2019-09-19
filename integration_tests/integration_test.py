import os
import subprocess
import time

from egcg_core import rest_communication, util, exceptions
from unittest.mock import Mock, patch

from egcg_core.app_logging import logging_default
from egcg_core.config import cfg
from egcg_core.integration_testing import ReportingAppIntegrationTest

from analysis_driver import client
from integration_tests import mocked_data


class IntegrationTest(ReportingAppIntegrationTest):
    patches = [
        patch('analysis_driver.client.load_config'),
        patch('egcg_core.clarity.get_plate_id_and_well', new=mocked_data.fake_get_plate_id_and_well),
        patch('egcg_core.clarity.get_sample_sex'),
        patch('egcg_core.clarity.get_sample_genotype', return_value=set()),
        patch('egcg_core.clarity.get_samples_arrived_with', return_value=set()),
        patch('egcg_core.clarity.get_samples_genotyped_with', return_value=set()),
        patch('egcg_core.clarity.get_samples_sequenced_with', return_value=set()),
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
    ]

    def __init__(self, *args):
        super().__init__(*args)
        self.run_id = self.cfg['input_data']['run_id']
        self.rapid_run_id = self.cfg['input_data']['rapid_run_id']
        self.aborted_run_id = self.cfg['input_data']['aborted_run_id']
        self.barcode = self.cfg['input_data']['barcode']
        self.project_id = self.cfg['input_data']['project_id']
        self.sample_id = self.cfg['input_data']['sample_id']
        self.library_id = self.cfg['input_data']['library_id']
        self.samples_for_project = self.cfg['input_data']['samples_for_project']
        self.run_dir = os.path.dirname(os.getcwd())  # we're inside the checked out project, not the top level

    def setUp(self):
        # Set the LIMS data yaml file before the super().setUp so it is picked up durring the setup
        self.lims_data_yaml = os.path.join(os.path.dirname(__file__), 'data_for_clarity_lims.yaml')

        super().setUp()
        cfg.load_config_file(os.getenv('ANALYSISDRIVERCONFIG'), env_var='ANALYSISDRIVERENV')

        run_elements = []
        for lane in range(1, 8):
            run_elements.append({
                'run_id': self.run_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
                'library_id': self.library_id, 'run_element_id': '%s_%s_%s' % (self.run_id, lane, self.barcode),
                'useable': 'yes', 'barcode': self.barcode, 'lane': lane, 'bases_r1': 70000000,
                'bases_r2': 70000000,'clean_bases_r1': 66000000, 'clean_bases_r2': 66000000,
                'q30_bases_r1': 64000000, 'q30_bases_r2': 64000000, 'clean_q30_bases_r1': 64000000,
                'clean_q30_bases_r2': 64000000,'clean_reads': 1
            })
        # Lane 8 has no data
        run_elements.append(
            {'run_id': self.run_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
             'library_id': self.library_id, 'run_element_id': '%s_%s_%s' % (self.run_id, 8, self.barcode),
             'useable': 'no', 'barcode': self.barcode, 'lane': 8, 'bases_r1': 0, 'bases_r2': 0,
             'clean_bases_r1': 0, 'clean_bases_r2': 0, 'q30_bases_r1': 0,
             'q30_bases_r2': 0, 'clean_q30_bases_r1': 0, 'clean_q30_bases_r2': 0,
             'clean_reads': 0}
        )
        for e in run_elements:
            rest_communication.post_entry('run_elements', e)

        rest_communication.post_entry(
            'samples',
            {'library_id': self.library_id, 'project_id': self.project_id, 'sample_id': self.sample_id,
             'run_elements': [e['run_element_id'] for e in run_elements], 'required_yield': 900000000,
             'required_coverage': 30}
        )

        rest_communication.post_entry('projects', {'project_id': self.project_id, 'samples': [self.sample_id]})

        rest_communication.post_entry('species', {'name': 'Homo sapiens', 'default_version': 'hg38'})
        rest_communication.post_entry('species', {'name': 'Canis lupus familiaris', 'default_version': 'CanFam3.1'})

        rest_communication.post_entry('genomes', {
            'assembly_name': 'hg38',
            'snpEff': 'GRCh38.86',
            'data_files': {
                'fasta': 'Homo_sapiens/hg38.fa',
                'variation': 'Homo_sapiens/hg38/dbsnp-147.vcf.gz',
                'vqsr': {
                    'hapmap': 'Homo_sapiens/hg38/gatk_bundle/hapmap_3.3.hg38.vcf.gz',
                    'omni': 'Homo_sapiens/hg38/gatk_bundle/1000G_omni2.5.hg38.vcf.gz',
                    'thousand_genomes': 'Homo_sapiens/hg38/gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz',
                    'dbsnp': 'Homo_sapiens/hg38/gatk_bundle/dbsnp_146.hg38.vcf.gz',
                    'mills': 'Homo_sapiens/hg38/gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
                }
            }
        })
        rest_communication.post_entry('genomes', {'assembly_name': 'CanFam3.1',
                                                  'data_files': {'fasta': 'Homo_sapiens/hg38.fa',
                                                                 'variation': 'Homo_sapiens/hg38/dbsnp-147.vcf.gz'}})
        rest_communication.post_entry('genomes', {'assembly_name': 'phix174',
                                                  'data_files': {
                                                      'fasta': 'PhiX/GCA_000819615.1_ViralProj14015_genomic.fna'}})

        self.dynamic_patches = []
        self._test_success = True

    def tearDown(self):
        super().tearDown()
        self._reset_logging()

    def setup_test(self, test_type, test_name, integration_section, species='Homo sapiens',
                   analysis_type='Variant Calling gatk'):
        cfg.content['jobs_dir'] = os.path.join(os.path.dirname(os.getcwd()), 'jobs', test_name)
        cfg.content[test_type]['output_dir'] = os.path.join(os.path.dirname(os.getcwd()), 'outputs', test_name)

        input_dir = cfg[test_type]['input_dir']
        if input_dir and os.path.isdir(input_dir):
            cfg.content[test_type]['input_dir'] = input_dir
        else:
            raise exceptions.EGCGError('Input dir %s does not exist' % input_dir)

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

        self._add_patches(
            patch('egcg_core.clarity.get_sample', new=_fake_get_sample),
            patch('egcg_core.clarity.get_species_from_sample', return_value=species)
        )

    @staticmethod
    def _reset_logging():
        for logger in logging_default.loggers.values():
            for handler in logger.handlers:
                logger.removeHandler(handler)

        logging_default.handlers = set()
        logging_default.loggers = {}

    def _add_patches(self, *patches):
        for p in patches:
            self.patches.append(p)
            p.start()

    def expect_equal(self, obs, exp, name):
        try:
            self.assertEqual(name, obs, exp)
        except AssertionError:
            self._test_success = False

    def expect_output_files(self, exp, base_dir):
        for k, v in exp.items():
            grepv = None
            if isinstance(v, dict):  # {'md5': 'an_md5sum', 'grepv': 'GATKCommandLine'}
                grepv = v.get('grepv')
                assert isinstance(grepv, list), 'grepv field for %s must be a list' % k
                v = v['md5']

            if base_dir:
                k = os.path.join(base_dir, k)

            if v is None:
                self.expect_equal(os.path.isfile(k), True, k)
            else:
                self.expect_equal(self._check_md5(k, grepv), v, k)

    def expect_qc_data(self, obs, exp):
        for e in sorted(exp):
            self.expect_equal(util.query_dict(obs, e), exp[e], e)

    def expect_stage_data(self, stage_names, **query_kw):
        """
        Take a list of expected stages, and compare with observed stages and exit statuses from the
        analysis_driver_stages endpoint.
        :param list stage_names: Expected stages and exit statuses. Each item can be a string stage name (in which case
                                 expected exit status will be 0), or a tuple of stage_name and exit_status.
        :param query_kw: Request parameters to pass to rest_communication.get_documents
        """
        stages = rest_communication.get_documents('analysis_driver_stages', **query_kw)
        obs = {s['stage_name']: s.get('exit_status') for s in stages}
        exp = dict([s if type(s) == tuple else (s, 0) for s in stage_names])
        self.expect_equal(obs, exp, 'stages')

    def _check_md5(self, fp, grepv_patterns=None):
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

        if grepv_patterns:
            for p in grepv_patterns:
                cmd += " | grep -v '{pattern}'".format(pattern=p)

        cmd += " | sed 's/{cwd}//g' | md5sum".format(cwd=self.run_dir.replace('/', '\/'))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        out, err = p.communicate()

        if err:
            raise ValueError(err)

        return out.decode('utf-8').rstrip('\n').split()[0]

    def test_demultiplexing(self):
        self.setup_test('run', 'test_demultiplexing', 'demultiplexing')
        self._add_patches(
            patch('analysis_driver.quality_control.interop_metrics.get_last_cycles_with_existing_bcls', return_value=310)
        )
        # Remove the Samplesheet that might have been generated before
        samplesheet = os.path.join(cfg['run']['input_dir'], self.run_id, 'SampleSheet_analysis_driver.csv')
        if os.path.isfile(samplesheet):
            os.remove(samplesheet)

        # Force the aborted run to be the first on in line
        rest_communication.post_entry('runs', {'run_id': self.run_id})
        payload = {
            'dataset_name': self.run_id,
            'proc_id': 'run_' + self.run_id + '_atime',
            'dataset_type': 'run',
            'status': 'force_ready'
        }
        rest_communication.post_or_patch('analysis_driver_procs', [payload], id_field='proc_id')
        # Ensure the next analysis_driver_procs won't be created at the same second.
        time.sleep(1)
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
        self.expect_stage_data([
            'setup', 'wellduplicates', 'bcl2fastq', 'phixdetection', 'fastqfilter', 'seqtkfqchk', 'md5sum', 'fastqc',
            'integritycheck', 'qcoutput1', 'dataoutput', 'cleanup', 'samtoolsdepthmulti', 'picardinsertsizemulti',
            'qcoutput2', 'runreview', 'picardmarkduplicatemulti', 'samtoolsstatsmulti', 'bwaalignmulti', 'waitforread2',
            'bcl2fastqpartialrun', 'picardgcbias'
        ])

        proc = rest_communication.get_document('analysis_driver_procs')
        self.expect_equal(
            rest_communication.get_document('runs', where={'run_id': self.run_id}).get('analysis_driver_procs'),
            [proc['proc_id']],
            'run proc registered'
        )
        self.expect_equal(
            proc['pipeline_used'],
            {'name': 'demultiplexing', 'toolset_version': 1, 'toolset_type': 'run_processing'},
            'pipeline used'
        )
        assert self._test_success

    def test_demultiplexing_aborted(self):
        self.setup_test('run', 'test_demultiplexing_aborted', 'demultiplexing')
        self._add_patches(
            patch('analysis_driver.quality_control.interop_metrics.get_last_cycles_with_existing_bcls', return_value=310)
        )
        # Force the aborted run to be the first on in line
        rest_communication.post_entry('runs', {'run_id': self.aborted_run_id})
        payload = {
            'dataset_name': self.aborted_run_id,
            'proc_id': 'run_' + self.aborted_run_id + '_atime',
            'dataset_type': 'run',
            'status': 'force_ready'
        }
        rest_communication.post_or_patch('analysis_driver_procs', [payload], id_field='proc_id')
        # Ensure the next analysis_driver_procs won't be created at the same second.
        time.sleep(1)

        exit_status = client.main(['--run'])

        self.assertEqual('exit status', exit_status, 0)
        ad_proc = rest_communication.get_document('analysis_driver_procs', where={'dataset_name': self.aborted_run_id}, sort='-_created')

        self.expect_equal(
            ad_proc['status'], 'aborted', 'pipeline status'
        )
        stages = [('setup', 9), ('waitforread2', 9)]
        self.expect_stage_data(stages, where={'analysis_driver_proc': ad_proc['proc_id']})
        assert self._test_success

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
                                'fixunmapped', 'blast', 'sexvalidation', 'vcfstats', 'samtoolsdepth',
                                'verifybamid', 'sampledataoutput', 'md5sum', 'cleanup', 'samplereview'])

        ad_proc = rest_communication.get_document('analysis_driver_procs')

        self.expect_equal(
            ad_proc['pipeline_used'],
            {'toolset_type': 'human_sample_processing', 'name': 'bcbio', 'toolset_version': 0},
            'pipeline used'
        )
        self.expect_equal(ad_proc['genome_used'], 'hg38', 'genome used')

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 8)],
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

        ad_proc = rest_communication.get_document('analysis_driver_procs')

        self.expect_equal(
            ad_proc['pipeline_used'],
            {'toolset_type': 'non_human_sample_processing', 'name': 'variant_calling', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(ad_proc['genome_used'], 'CanFam3.1', 'genome used')

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 8)],
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

        ad_proc = rest_communication.get_document('analysis_driver_procs')

        self.expect_equal(
            ad_proc['pipeline_used'],
            {'toolset_type': 'non_human_sample_processing', 'name': 'qc', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(ad_proc['genome_used'], 'CanFam3.1', 'genome used')

        self.expect_equal(
            rest_communication.get_document('analysis_driver_procs')['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 8)],
            'data source'
        )

    def test_qc(self):
        self.setup_test('sample', 'test_qc', 'qc', 'Canis lupus familiaris', 'QC GATK3')

        self._run_qc_test()
        self.expect_output_files(
            self.cfg['qc']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        assert self._test_success

    def test_resume(self):
        self.setup_test('sample', 'test_resume', 'qc', 'Canis lupus familiaris', 'QC GATK3')
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

        self.expect_stage_data(['genotypegvcfs', 'relatedness', 'peddy', 'parserelatedness', 'output', 'cleanup'])
        ad_procs = rest_communication.get_document('analysis_driver_procs', where={'dataset_name': self.project_id})
        self.expect_equal(
            ad_procs['pipeline_used'],
            {'toolset_type': 'project_processing', 'name': 'project', 'toolset_version': 0},
            'pipeline used'
        )

        assert self._test_success

    def test_rapid_analysis(self):
        self._add_patches(
            patch('egcg_core.clarity.get_run', return_value=mocked_data.mocked_rapid_run),
            patch('analysis_driver.quality_control.interop_metrics.get_last_cycles_with_existing_bcls', return_value=302)
        )

        sample_ids = ['non_pooling_sample_%s' % i for i in range(1, 9)]
        for s in sample_ids:
            rest_communication.post_entry(
                'samples', {'project_id': 'rapid_project', 'sample_id': s, 'user_sample_id': 'uid_' + s}
            )

        rest_communication.patch_entry(
            'projects',
            {'samples': sample_ids},
            'project_id',
            self.project_id,
            update_lists=['samples']
        )

        self.setup_test('run', 'test_rapid', 'rapid')
        # Remove the Samplesheet that might have been generated before
        samplesheet = os.path.join(cfg['run']['input_dir'], self.rapid_run_id, 'SampleSheet_analysis_driver.csv')
        if os.path.isfile(samplesheet):
            os.remove(samplesheet)

        cfg.content['delivery'] = {'dest': os.path.join(os.path.dirname(os.getcwd()), 'delivered_outputs', 'test_rapid')}
        os.makedirs(cfg['delivery']['dest'])
        cfg.content['sample']['output_dir'] = os.path.join(os.path.dirname(os.getcwd()), 'outputs', 'test_rapid')

        exit_status = client.main(['--run'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_output_files(
            self.cfg['rapid']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id)
        )
        self.expect_output_files(
            self.cfg['rapid']['delivered_files'],
            base_dir=util.find_file(cfg['delivery']['dest'], self.project_id, '*')
        )
        self.expect_qc_data(
            rest_communication.get_document(
                'samples',
                where={'sample_id': 'non_pooling_sample_2'}
            ),
            self.cfg['rapid']['qc']
        )
        self.expect_stage_data(['setup', 'wellduplicates', 'bcl2fastq', 'phixdetection', 'fastqfilter', 'seqtkfqchk',
                                'md5sum', 'fastqc', 'integritycheck', 'qcoutput1', 'dataoutput', 'cleanup',
                                'samtoolsdepthmulti', 'picardinsertsizemulti', 'qcoutput2', 'runreview',
                                'picardmarkduplicatemulti', 'samtoolsstatsmulti', 'bwaalignmulti', 'waitforread2',
                                'bcl2fastqpartialrun', 'picardgcbias', 'dragen', 'dragenmetrics', 'dragenoutput'])

        proc = rest_communication.get_document('analysis_driver_procs')
        self.expect_equal(
            rest_communication.get_document('runs', where={'run_id': self.rapid_run_id}).get('analysis_driver_procs'),
            [proc['proc_id']],
            'run proc registered'
        )
        self.expect_equal(
            proc['pipeline_used'],
            {'name': 'demultiplexing', 'toolset_version': 1, 'toolset_type': 'run_processing'},
            'pipeline used'
        )
        assert self._test_success

    def test_gatk4_qc(self):
        self.setup_test('sample', 'test_gatk4_qc', 'gatk4_qc', 'Canis lupus familiaris', 'Not Variant Calling')

        run_elements = []
        # remove  Yield from the run elements and add coverage so it starts because it passes the coverage threshold.
        for lane in range(1, 8):
            run_elements.append({
                'run_element_id': '%s_%s_%s' % (self.run_id, lane, self.barcode),
                'bases_r1': 1, 'bases_r2': 1, 'clean_bases_r1': 1, 'clean_bases_r2': 1,
                'q30_bases_r1': 1, 'q30_bases_r2': 1, 'clean_q30_bases_r1': 1, 'clean_q30_bases_r2': 1,
                'coverage': {'mean': 5}  # 7 * 5 = 35X coverage
            })

        exit_status = client.main(['--sample'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_qc_data(
            rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
            self.cfg['gatk4_qc']['qc']
        )

        self.expect_output_files(
            self.cfg['gatk4_qc']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )
        self.expect_stage_data(['mergefastqs', 'fastqscreen', 'blast', 'fastqindex', 'splitbwa', 'mergebamanddup',
                                'samtoolsstats', 'samtoolsdepth', 'splithaplotypecaller', 'gathervcf', 'selectsnps',
                                'selectindels', 'snpsfiltration', 'indelsfiltration', 'mergevariants', 'vcfstats',
                                'sampledataoutput', 'md5sum', 'samplereview', 'cleanup'])

        ad_proc = rest_communication.get_document('analysis_driver_procs')

        self.expect_equal(
            ad_proc['pipeline_used'],
            {'toolset_type': 'gatk4_sample_processing', 'name': 'qc_gatk4', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(ad_proc['genome_used'], 'CanFam3.1', 'genome used')

        self.expect_equal(
            ad_proc['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 8)],
            'data source'
        )

        assert self._test_success

    def test_gatk4_var_calling(self):
        self.setup_test(test_type='sample', test_name='test_gatk4_var_calling',
                        integration_section='gatk4_var_calling', species='Canis lupus familiaris',
                        analysis_type='Variant Calling gatk4')
        exit_status = client.main(['--sample'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_qc_data(
            rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
            self.cfg['gatk4_var_calling']['qc']
        )

        self.expect_output_files(
            self.cfg['gatk4_var_calling']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        self.expect_stage_data([
            'gathervcfvc', 'scatterbaserecalibrator', 'samplereview', 'vcfstats', 'indelsfiltration', 'md5sum',
            'selectindels', 'cleanup', 'splitgenotypegvcfs', 'fastqindex', 'gatherrecalbam',
            'merge_variants_hard_filter', 'verifybamid', 'samtoolsstats', 'gatherbqsrreport',
            'gathergvcf', 'splithaplotypecallervc', 'fastqscreen', 'selectsnps', 'snpsfiltration', 'splitbwa',
            'blast', 'scatterapplybqsr', 'mergefastqs', 'samtoolsdepth', 'sampledataoutput', 'mergebamanddup'
        ])

        ad_proc = rest_communication.get_document('analysis_driver_procs')

        self.expect_equal(ad_proc['genome_used'], 'CanFam3.1', 'genome used')

        self.expect_equal(
            ad_proc['pipeline_used'],
            {'toolset_type': 'gatk4_sample_processing', 'name': 'variant_calling_gatk4', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(
            ad_proc['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 8)],
            'data source'
        )
        assert self._test_success

    def test_gatk4_var_calling_human(self):
        self.setup_test(test_type='sample', test_name='test_gatk4_var_calling_human',
                        integration_section='gatk4_var_calling_human',
                        species='Homo sapiens', analysis_type='Variant Calling gatk4')
        exit_status = client.main(['--sample'])
        self.assertEqual('exit status', exit_status, 0)

        self.expect_qc_data(
            rest_communication.get_document('samples', where={'sample_id': self.sample_id}),
            self.cfg['gatk4_var_calling_human']['qc']
        )

        self.expect_output_files(
            self.cfg['gatk4_var_calling_human']['files'],
            base_dir=os.path.join(cfg['sample']['output_dir'], self.project_id, self.sample_id)
        )

        self.expect_stage_data([
            'gathervcfvc', 'mergebamanddup', 'splitgenotypegvcfs', 'selectsnps', 'mergefastqs', 'cleanup',
            'splithaplotypecallervc', 'variantannotation', 'genotypevalidation', 'gatherbqsrreport',
            'selectindels', 'verifybamid', 'gathergvcf', 'sexvalidation', 'fastqscreen', 'sampledataoutput',
            'gatherrecalbam', 'indelsfiltration', 'samtoolsdepth', 'scatterapplybqsr', 'samplereview', 'fastqindex',
            'scatterbaserecalibrator', 'merge_variants_hard_filter', 'blast', 'splitbwa', 'vcfstats', 'md5sum',
            'samtoolsstats', 'snpsfiltration'
        ])

        ad_proc = rest_communication.get_document('analysis_driver_procs')

        self.expect_equal(ad_proc['genome_used'], 'hg38', 'genome used')

        self.expect_equal(
            ad_proc['pipeline_used'],
            {'toolset_type': 'gatk4_sample_processing', 'name': 'human_variant_calling_gatk4', 'toolset_version': 0},
            'pipeline used'
        )

        self.expect_equal(
            ad_proc['data_source'],
            ['_'.join([self.run_id, str(i), self.barcode]) for i in range(1, 8)],
            'data source'
        )

        assert self._test_success
