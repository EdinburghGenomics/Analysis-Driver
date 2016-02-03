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
        reference = self.validator.validation_cfg['reference']
        for sample_id in ['10015AT0001', '10015AT0002']:
            cmd = self.validator._bwa_aln(
                self.mapping[sample_id],
                sample_id,
                sample_id + '.bam',
                cfg.query('genotype-validation', 'reference')
            )
            expected = (
                cfg['tools']['bwa'],
                'sampe -r \'@RG\\tID:1\\tSM:%s\'' % sample_id,
                reference,
                '<(' + cfg['tools']['bwa'], 'aln', reference, self.mapping[sample_id][0] + ')',
                '<(' + cfg['tools']['bwa'], 'aln', reference, self.mapping[sample_id][1] + ')',
                self.mapping[sample_id][0],
                self.mapping[sample_id][1],
                '|',
                cfg['tools']['samblaster'], '--removeDups',
                '|',
                cfg['tools']['samtools'], 'view -F 4 -Sb -',
                '|',
                cfg['tools']['sambamba'], 'sort -t 16 -o ', sample_id + '.bam', '/dev/stdin'
            )
            self.compare_lists(cmd, ' '.join(expected))
