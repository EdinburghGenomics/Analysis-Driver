from tests.test_analysisdriver import TestAnalysisDriver
import os.path
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import GenotypeValidation
__author__ = 'tcezard'


class TestGenotypeValidation(TestAnalysisDriver):
    @property
    def mapping(self):
        return {
            '10015AT0001': [
                os.path.join(self.fastq_path, '10015AT', '10015AT0001', x) for x in [
                    '10015AT0001_S6_L004_R1_001.fastq.gz', '10015ATA0001_merged_R2.fastq.gz'
                ]
            ],
            '10015AT0002': [
                os.path.join(self.fastq_path, '10015AT', '10015AT0002', x) for x in [
                    '10015AT0002_merged_R1.fastq.gz', '10015ATA0002_merged_R2.fastq.gz'
                ]
            ]
        }

    def setUp(self):
        self.validator = GenotypeValidation(self.mapping, 'test_run')

    def test_bwa_aln(self):
        vc = self.validator.validation_cfg
        for sample_id in ['10015AT0001', '10015AT0002']:
            cmd = self.validator._bwa_aln(
                self.mapping[sample_id],
                sample_id,
                sample_id + '.bam',
                cfg.query('genotype-validation', 'reference')
            )
            expected = (
                vc['bwa'],
                'sampe -r \'@RG\\tID:1\\tSM:%s\'' % sample_id,
                vc['reference'],
                '<(' + vc['bwa'], 'aln', vc['reference'], self.mapping[sample_id][0] + ')',
                '<(' + vc['bwa'], 'aln', vc['reference'], self.mapping[sample_id][1] + ')',
                self.mapping[sample_id][0],
                self.mapping[sample_id][1],
                '|',
                vc['samblaster'], '--removeDups',
                '|',
                vc['samtools'], 'view -F 4 -Sb -',
                '|',
                vc['sambamba'], 'sort -t 16 -o ', sample_id + '.bam', '/dev/stdin'
            )
            assert cmd == ' '.join(expected)
