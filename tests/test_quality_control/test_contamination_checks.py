from os.path import join
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import FastqScreen, VerifyBamID, Blast
from unittest.mock import Mock, patch

ppath = 'analysis_driver.quality_control.contamination_checks.'


class TestContaminationCheck(QCTester):
    def setUp(self):
        super().setUp()
        self.c = FastqScreen(dataset=self.dataset, fq_pattern='r?.fastq.gz')

    def test_fastqscreen_command(self):
        with patch(ppath + 'util.find_files', return_value=['r1.fastq.gz', 'r2.fastq.gz']):
            assert self.c._fastqscreen_command() == (
                'path/to/fastqscreen --aligner bowtie2 r1.fastq.gz r2.fastq.gz '
                '--conf path/to/fastqscreen/conf --force'
            )

    def test_get_expected_outfiles(self):
        assert self.c.fastqscreen_expected_outfiles == 'r?_screen.txt'

    @patch('egcg_core.executor.execute')
    def test_run_fastqscreen(self, mocked_execute):
        with patch(ppath + 'FastqScreen._fastqscreen_command', return_value='a_cmd'):
            self.c._run()

        mocked_execute.assert_called_once_with(
            'a_cmd', working_dir='tests/assets/jobs/test_sample', mem=10, cpus=2, job_name='fastqscreen'
        )


class TestVerifyBamID(QCTester):
    def setUp(self):
        super().setUp()
        self.vbi = VerifyBamID(dataset=self.dataset, bam_file='test_bam_file.bam')

    @patch('egcg_core.executor.execute')
    def test_contamination_check(self, mocked_execute):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            self.vbi._run()

        mocked_execute.assert_any_call(
            'path/to/samtools view -b test_bam_file.bam chr22 > tests/assets/jobs/test_sample/test_sample_chr22.bam',
            mem=2,
            job_name='filter_bam22',
            working_dir='tests/assets/jobs/test_sample',
            cpus=1,
            log_commands=False
        )
        mocked_execute.assert_any_call(
            'path/to/samtools index tests/assets/jobs/test_sample/test_sample_chr22.bam',
            mem=2,
            working_dir='tests/assets/jobs/test_sample',
            cpus=1,
            job_name='index_bam22'
        )
        mocked_execute.assert_any_call(
            'path/to/verifybamid --bam tests/assets/jobs/test_sample/test_sample_chr22.bam --vcf path/to/population_vcf --out tests/assets/jobs/test_sample/test_sample-chr22-vbi',
            working_dir='tests/assets/jobs/test_sample',
            cpus=1,
            job_name='verify_bam_id',
            mem=4
        )


tax_dict = {40674: 'Mammalia', 9443: 'Primates', 9604: 'Hominidae', 9605: 'Homo', 9606: 'Homo sapiens',
            33208: 'Metazoa', 7711: 'Chordata', 2759: 'Eukaryota', 99802: 'Spirometra erinaceieuropaei',
            33154: 'Opisthokonta', 28843: 'Diphyllobothriidae', 6157: 'Platyhelminthes',
            131567: 'cellular organisms', 6072: 'Eumetazoa', 46580: 'Spirometra', 6199: 'Cestoda',
            6200: 'Eucestoda', 33213: 'Bilateria', 1224679: 'Diphyllobothriidea'}

homo_sapiens_tree = {
    'Eukaryota': {
        'reads': 197,
        'Metazoa': {
            'reads': 197,
            'Chordata': {
                'reads': 197,
                'Mammalia': {
                    'reads': 197,
                    'Primates': {
                        'reads': 197,
                        'Hominidae': {
                            'reads': 197,
                            'Homo': {
                                'reads': 197,
                                'Homo sapiens': {
                                    'reads': 197
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


def mocked_txid_transl(taxids):
    return {taxid: tax_dict.get(taxid, 'Unavailable') for taxid in taxids}


class TestBlast(QCTester):
    def setUp(self):
        super().setUp()
        self.blast = Blast(dataset=self.run_dataset, fastq_file='fastqFile1.fastq')
        self.blast._ncbi = Mock(get_taxid_translator=mocked_txid_transl)

    def test_sample_fastq_command(self):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            assert self.blast.sample_fastq_command() == (
                'set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | '
                'path/to/seqtk seq -a > tests/assets/jobs/test_run/fastqFile1_sample3000.fasta'
            )

    def test_fasta_blast_command(self):
        assert self.blast.fasta_blast_command() == (
            'export PATH=$PATH:/path/to/db/dir; path/to/blastn '
            '-query tests/assets/jobs/test_run/fastqFile1_sample3000.fasta -db path/to/db/dir/nt '
            '-out tests/assets/jobs/test_run/fastqFile1_sample3000_blastn -num_threads 12 -max_target_seqs 1 -max_hsps 1 '
            "-outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"
        )

    def test_get_taxids(self):
        blast_file = join(self.assets_path, 'blast_outfile')
        taxids = self.blast.get_taxids(blast_file)
        assert taxids == {'9598': 2, '9606': 7, '99802': 2}

    def test_get_ranks(self):
        expected_res = {9604: 'family', 33208: 'kingdom', 89593: 'subphylum'}
        self.blast._ncbi.get_rank = Mock(return_value=expected_res)
        self.blast._ncbi.get_lineage = Mock(return_value=[9604, 33208, 89593])

        assert self.blast.get_ranks('9960') == expected_res
        self.blast.ncbi.get_lineage.assert_called_once_with(9960)
        self.blast.ncbi.get_rank.assert_called_once_with([9604, 33208, 89593])

        self.blast._ncbi.get_lineage.side_effect = ValueError
        assert self.blast.get_ranks('9960') == {0: 'rank unavailable'}

    @patch('analysis_driver.quality_control.Blast.get_ranks')
    def test_get_all_taxa_identified1(self, mocked_rank):
        mocked_rank.return_value = {
            9604: 'family', 9605: 'genus', 9606: 'species', 2759: 'superkingdom', 376913: 'suborder',
            314146: 'superorder', 7711: 'phylum', 40674: 'class', 9443: 'order', 207598: 'subfamily',
            314293: 'infraorder', 9526: 'parvorder', 314295: 'superfamily', 33208: 'kingdom',
            89593: 'subphylum'
        }
        observed = self.blast.get_all_taxa_identified(taxon_dict={}, taxon='9606', taxids={'9606': 197})
        assert observed == homo_sapiens_tree

    @patch('analysis_driver.quality_control.Blast.get_ranks')
    def test_get_all_taxa_identified2(self, mocked_rank):
        mocked_rank.return_value = {
            99802: 'species', 1: 'no rank', 33154: 'no rank', 2759: 'superkingdom', 28843: 'family',
            6157: 'phylum', 131567: 'no rank', 6072: 'no rank', 46580: 'genus', 6199: 'class',
            33208: 'kingdom', 6200: 'subclass', 33213: 'no rank', 1224679: 'order'
        }

        test_spirometra_and_human = self.blast.get_all_taxa_identified(
            taxon_dict=homo_sapiens_tree, taxon='99802', taxids={'99802': 2}
        )

        expected = homo_sapiens_tree.copy()
        expected['Eukaryota']['reads'] = 199
        expected['Eukaryota']['Metazoa']['reads'] = 199
        expected['Eukaryota']['Metazoa']['Platyhelminthes'] = {
            'reads': 2,
            'Cestoda': {
                'reads': 2,
                'Diphyllobothriidea': {
                    'reads': 2,
                    'Diphyllobothriidae': {
                        'reads': 2,
                        'Spirometra': {
                            'reads': 2,
                            'Spirometra erinaceieuropaei': {'reads': 2}
                        }
                    }
                }
            }
        }
        assert test_spirometra_and_human == expected

    @patch('egcg_core.executor.execute')
    def test_run_sample_fastq(self, mocked_execute):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            self.blast.run_sample_fastq()
        mocked_execute.assert_called_with(
            'set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | path/to/seqtk seq -a > tests/assets/jobs/test_run/fastqFile1_sample3000.fasta',
            working_dir='tests/assets/jobs/test_run',
            mem=10,
            cpus=2,
            job_name='sample_fastq'
        )

    @patch('egcg_core.executor.execute')
    def test_run_blast(self, mocked_execute):
        with patch(ppath + 'util.find_file', new=self.fake_find_file):
            self.blast.run_blast()
        mocked_execute.assert_called_with(
            ('export PATH=$PATH:/path/to/db/dir; path/to/blastn '
             '-query tests/assets/jobs/test_run/fastqFile1_sample3000.fasta -db path/to/db/dir/nt '
             '-out tests/assets/jobs/test_run/fastqFile1_sample3000_blastn -num_threads 12 -max_target_seqs 1 -max_hsps '
             "1 -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"),
            working_dir='tests/assets/jobs/test_run',
            mem=20,
            cpus=12,
            job_name='contamination_blast'
        )
