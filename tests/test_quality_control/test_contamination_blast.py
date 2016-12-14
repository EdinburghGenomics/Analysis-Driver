import os
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationBlast
from unittest.mock import patch, Mock

tax_dict = {
    40674: 'Mammalia',
    9443: 'Primates',
    9604: 'Hominidae',
    9605: 'Homo',
    9606: 'Homo sapiens',
    33208: 'Metazoa',
    7711: 'Chordata',
    2759: 'Eukaryota',
    99802: 'Spirometra erinaceieuropaei',
    33154: 'Opisthokonta',
    28843: 'Diphyllobothriidae',
    6157: 'Platyhelminthes',
    131567: 'cellular organisms',
    6072: 'Eumetazoa',
    46580: 'Spirometra',
    6199: 'Cestoda',
    6200: 'Eucestoda',
    33213: 'Bilateria',
    1224679: 'Diphyllobothriidea'
}

homo_sapiens_tree = {
    'Eukaryota': {
        'Metazoa': {
            'Chordata': {
                'reads': 197,
                'Mammalia': {
                    'Primates': {
                        'reads': 197,
                        'Hominidae': {
                            'reads': 197,
                            'Homo': {
                                'Homo sapiens': {
                                    'reads': 197
                                },
                                'reads': 197
                            }
                        }
                    },
                    'reads': 197
                }
            },
            'reads': 197
        },
        'reads': 197
    }
}


def mocked_txid_transl(taxids):
    return {taxid: tax_dict.get(taxid, 'Unavailable') for taxid in taxids}


class TestContaminationBlast(QCTester):
    def setUp(self):
        super().setUp()
        self.fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        self.fastq_file = 'fastqFile1.fastq'
        self.fasta_file = 'fastaFile1.fasta'
        self.working_dir = 'test_run'
        self.blast = ContaminationBlast(self.run_dataset, self.working_dir, self.fastq_files)
        self.blast._ncbi = Mock(get_taxid_translator=mocked_txid_transl)

    def test_sample_fastq_command(self):
        command, outfile = self.blast.sample_fastq_command(self.fastq_file, nb_reads=3000)
        assert command == ('set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | '
                           'path/to/seqtk seq -a > test_run/fastqFile1_sample3000.fasta')
        assert outfile == 'test_run/fastqFile1_sample3000.fasta'

    def test_fasta_blast_command(self):
        command, outfile = self.blast.fasta_blast_command(self.fasta_file)
        assert command == (
            'export PATH=$PATH:/path/to/db/dir; path/to/blastn -query fastaFile1.fasta -db path/to/db/dir/nt '
            '-out test_run/fastaFile1_blastn -num_threads 12 -max_target_seqs 1 -max_hsps 1 '
            "-outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"
        )
        assert outfile == 'test_run/fastaFile1_blastn'

    def test_get_taxids(self):
        blast_file = os.path.join(self.assets_path, 'blast_outfile')
        taxids = self.blast.get_taxids(blast_file)
        assert taxids == {'9598': 2, '9606': 7, '99802': 2}

    def test_get_ranks(self):
        expected_res = {9604: 'family', 33208: 'kingdom', 89593: 'subphylum'}
        self.blast._ncbi.get_rank = Mock(return_value=expected_res)
        self.blast._ncbi.get_lineage = Mock(return_value=[9604, 33208, 89593])

        assert self.blast.get_ranks('9960') == expected_res
        self.blast.ncbi.get_lineage.assert_called_once_with(9960)
        self.blast.ncbi.get_rank.assert_called_once_with([9604, 33208, 89593])

        with patch('analysis_driver.quality_control.ContaminationBlast.ncbi', new_callable=PropertyMock) as ncbi:
            e = ValueError
            ncbi().get_lineage.side_effect = e
            assert self.contamination_blast.get_ranks('9960') == {0: 'rank unavailable'}

    @patch('analysis_driver.quality_control.ContaminationBlast.get_ranks')
    def test_get_all_taxa_identified1(self, mocked_rank):
        mocked_rank.return_value = {
            9604: 'family',
            9605: 'genus',
            9606: 'species',
            2759: 'superkingdom',
            376913: 'suborder',
            314146: 'superorder',
            7711: 'phylum',
            40674: 'class',
            9443: 'order',
            207598: 'subfamily',
            314293: 'infraorder',
            9526: 'parvorder',
            314295: 'superfamily',
            33208: 'kingdom',
            89593: 'subphylum'
        }
        observed = self.blast.get_all_taxa_identified(taxon_dict={}, taxon='9606', taxids={'9606': 197})
        assert observed == homo_sapiens_tree

    @patch('analysis_driver.quality_control.ContaminationBlast.get_ranks')
    def test_get_all_taxa_identified2(self, mocked_rank):
        mocked_rank.return_value = {
            99802: 'species',
            1: 'no rank',
            33154: 'no rank',
            2759: 'superkingdom',
            28843: 'family',
            6157: 'phylum',
            131567: 'no rank',
            6072: 'no rank',
            46580: 'genus',
            6199: 'class',
            33208: 'kingdom',
            6200: 'subclass',
            33213: 'no rank',
            1224679: 'order'
        }

        test_spirometra_and_human = self.blast.get_all_taxa_identified(
            taxon_dict=homo_sapiens_tree, taxon='99802', taxids={'99802': 2}
        )

        expected = dict(homo_sapiens_tree)
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

    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_sample_fastq(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        sample_fastq_outfile = self.blast.run_sample_fastq(self.fastq_file, nb_reads=3000)
        assert sample_fastq_outfile == 'test_run/fastqFile1_sample3000.fasta'
        mocked_execute.assert_called_once_with(
            'set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | path/to/seqtk seq -a > test_run/fastqFile1_sample3000.fasta',
            working_dir='test_run',
            mem=10,
            cpus=2,
            job_name='sample_fastq'
        )

    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_blast(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        blast_outfile = self.blast.run_blast(self.fasta_file)
        assert blast_outfile == 'test_run/fastaFile1_blastn'
        mocked_execute.assert_called_once_with(
            ('export PATH=$PATH:/path/to/db/dir; path/to/blastn -query fastaFile1.fasta -db '
             'path/to/db/dir/nt -out test_run/fastaFile1_blastn -num_threads 12 -max_target_seqs 1 -max_hsps '
             "1 -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"),
            working_dir='test_run',
            mem=20,
            cpus=12,
            job_name='contamination_blast'
        )
