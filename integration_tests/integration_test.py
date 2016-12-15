import os
import sys
import gzip
import hashlib
import pytest
from io import StringIO
from datetime import datetime
from contextlib import redirect_stdout, contextmanager
from egcg_core.config import cfg, Configuration
from egcg_core.app_logging import logging_default
from egcg_core import rest_communication, notifications, util, executor
from unittest.mock import Mock, patch
from analysis_driver import client

cfg.load_config_file(os.getenv('ANALYSISDRIVERCONFIG'), env_var='ANALYSISDRIVERENV')
integration_cfg = Configuration(os.getenv('INTEGRATIONCONFIG'))
app_logger = logging_default.get_logger(__name__)
comm = rest_communication.Communicator((cfg['rest_api']['username'], cfg['rest_api']['password']), cfg['rest_api']['url'])

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
def patch_pipeline(species='Homo sapiens', analysis_type='Variant Calling'):
    patches = []

    def _patch(ppath, **kwargs):
        p = patch('analysis_driver.' + ppath, **kwargs)
        p.start()
        patches.append(p)

    def _fake_get_sample(sample_name):
        return Mock(name=sample_name, udf={'Coverage': 1337, 'Analysis Type': analysis_type})

    _patch('driver.clarity.get_species_from_sample', return_value=species)
    _patch('driver.clarity.get_sample', new=_fake_get_sample)
    _patch('report_generation.report_crawlers.clarity.get_species_from_sample', return_value=species)
    _patch('report_generation.report_crawlers.clarity.get_sample', new=_fake_get_sample)
    _patch('dataset.get_expected_yield_for_sample', return_value=0.9)
    _patch('report_generation.report_crawlers.clarity.get_expected_yield_for_sample', return_value=0.9)
    _patch('dataset_scanner.get_list_of_samples', new=_fake_get_list_of_samples)
    _patch('driver.clarity.get_run', return_value=Mock(udf={'Run Status': 'RunCompleted'}))
    _patch('driver.clarity.find_project_name_from_sample', return_value='a_project')
    _patch('quality_control.genotype_validation.clarity.find_project_name_from_sample', return_value='a_project')
    _patch('driver.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('report_generation.report_crawlers.clarity.get_user_sample_name', new=_fake_get_user_sample_id)
    _patch('report_generation.report_crawlers.clarity.get_plate_id_and_well', new=_fake_get_plate_id_and_well)
    _patch('report_generation.report_crawlers.clarity.get_sample_gender')
    _patch('quality_control.genotype_validation.clarity.get_samples_arrived_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_samples_genotyped_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_samples_sequenced_with', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_sample_names_from_project', return_value=set())
    _patch('quality_control.genotype_validation.clarity.get_sample_genotype', return_value=set())
    _patch('quality_control.lane_duplicates.WellDuplicates._well_duplicates', new=_fake_welldups)
    _patch('driver.time.sleep')

    yield

    for p in patches:
        p.stop()


def test_demultiplexing():
    with patch_pipeline():
        sys.argv = [entry_point, '--run']
        exit_status = client.main()
        assert exit_status == 0

        comm.patch_entries('run_elements', {'useable': 'yes'}, where={'sample_id': '10015AT0004'})
        assert sorted(comm.get_document('projects')['samples']) == ['10015AT000' + str(i) for i in (1, 2, 3, 4, 6, 7, 8, 9)]
        assert len(comm.get_document('samples', where={'sample_id': '10015AT0004'})['run_elements']) == 8
        output_dir = os.path.join(cfg['run']['output_dir'], '150723_E00306_0025_BHCHK3CCXX')
        output_fastqs = util.find_files(output_dir, '*.fastq.gz') + util.find_files(output_dir, '10015AT', '*', '*.fastq.gz')
        assert len(output_fastqs) == 126  # 14 undetermined + 112 sample

        for fq in output_fastqs:
            with gzip.open(fq, 'r') as f:
                m = hashlib.md5()
                m.update(f.read())
                assert m.hexdigest() == integration_cfg['demultiplexing']['md5s'][os.path.basename(fq)]


def test_bcbio():
    with patch_pipeline():
        sys.argv = [entry_point, '--sample']
        exit_status = client.main()
        assert exit_status == 0

        e = rest_communication.get_document('samples', where={'sample_id': '10015AT0004'})
        for obs, exp in ((e['bam_file_reads'], integration_cfg['bcbio']['qc']['bam_file_reads']),
                         (e['coverage']['genome_size'], integration_cfg['bcbio']['qc']['genome_size']),
                         (e['duplicate_reads'], integration_cfg['bcbio']['qc']['duplicate_reads']),
                         (e['called_gender'], integration_cfg['bcbio']['qc']['called_gender'])):
            assert obs == exp

        output_dir = os.path.join(cfg['sample']['output_dir'], 'a_project', '10015AT0004')

        for outfile, md5file in (('uid_10015AT0004.txt', 'samtools_stats.txt.md5'),
                                 ('uid_10015AT0004.vcf.stats', 'uid_10015AT0004.vcf.stats.md5')):
            with open(os.path.join(output_dir, md5file), 'r') as f:
                assert f.readline().split(' ')[0] == integration_cfg['bcbio']['md5s'][outfile]


# TODO: var_calling_pipeline, qc_pipeline


def main():
    start_time = _now()
    s = StringIO()
    with redirect_stdout(s):
        pytest.main([__file__])
    end_time = _now()

    test_output = util.str_join(
        'Pipeline test finished. ',
        'Start time: %s, finish time: %s. '
        'Pytest output below:\n' % (start_time, end_time),
        s.getvalue()
    )
    e = notifications.EmailNotification(
        'Analysis Driver integration test',
        **integration_cfg['notification']
    )
    e.notify(test_output)


if __name__ == '__main__':
    main()
