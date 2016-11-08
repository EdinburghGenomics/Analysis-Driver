import os
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationBlast
from unittest.mock import patch
from collections import Counter

class TestContaminationBlast(QCTester):
    def setUp(self):
        super().setUp()
        self.fastq_files = ['fastqFile1.fastq', 'fastqFile2.fastq']
        self.fastq_file = 'fastqFile1.fastq'
        self.fasta_file = 'fastaFile1.fasta'
        self.working_dir = 'test_run'
        self.contamination_blast = ContaminationBlast(self.run_dataset, self.working_dir, self.fastq_files)

    def test_sample_fastq_command(self):
        command, outfile = self.contamination_blast.sample_fastq_command(self.fastq_file)
        assert command == 'set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | path/to/seqtk seq -a > test_run/fastqFile1_sample3000.fasta'
        assert outfile == 'test_run/fastqFile1_sample3000.fasta'

    def test_fasta_blast_command(self):
        command, outfile = self.contamination_blast.fasta_blast_command(self.fasta_file)
        assert command == "export PATH=$PATH:/path/to/db/dir; path/to/blastn -query fastaFile1.fasta -db path/to/db/dir/nt -out test_run/fastaFile1_blastn -num_threads 12 -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"
        assert outfile =='test_run/fastaFile1_blastn'

    def test_get_taxids(self):
        blast_file = os.path.join(self.assets_path, 'blast_outfile')
        taxids = self.contamination_blast.get_taxids(blast_file)
        assert taxids == {'9598': 2, '9606': 197, '99802': 2}

    @patch('analysis_driver.quality_control.ContaminationBlast.get_rank')
    @patch('analysis_driver.quality_control.ContaminationBlast.ncbi')
    @patch('analysis_driver.quality_control.ContaminationBlast.taxid_translator')
    def test_get_all_taxa_identified1(self, mocked_taxid_translator, mocked_ncbi, mocked_rank):
        mocked_rank.return_value = {1: 'no rank',
                                    33154: 'no rank',
                                    9347: 'no rank',
                                    9604: 'family',
                                    9605: 'genus',
                                    9606: 'species',
                                    2759: 'superkingdom',
                                    32523: 'no rank',
                                    32524: 'no rank',
                                    32525: 'no rank',
                                    376913: 'suborder',
                                    6072: 'no rank',
                                    1437010: 'no rank',
                                    117571: 'no rank',
                                    314146: 'superorder',
                                    117570: 'no rank',
                                    7711: 'phylum',
                                    7776: 'no rank',
                                    40674: 'class',
                                    9443: 'order',
                                    33511: 'no rank',
                                    207598: 'subfamily',
                                    131567: 'no rank',
                                    1338369: 'no rank',
                                    314293: 'infraorder',
                                    9526: 'parvorder',
                                    314295: 'superfamily',
                                    33208: 'kingdom',
                                    89593: 'subphylum',
                                    8287: 'no rank',
                                    33213: 'no rank',
                                    7742: 'no rank'}
        mocked_ncbi.return_value = None
        mocked_taxid_translator.return_value = ('Metazoa', 'Chordata', 'Mammalia', 'Primates', 'Hominidae', 'Homo', 'Homo sapiens')
        taxon_dict = {}
        taxon = '9606'
        taxids = Counter({'9606': 197, '99802': 2, '9598': 2})
        test_human_only = self.contamination_blast.get_all_taxa_identified(taxon_dict, taxon, taxids)
        assert test_human_only == {'Metazoa': {'Chordata': {'reads': 197, 'Mammalia': {'Primates': {'reads': 197, 'Hominidae': {'reads': 197, 'Homo': {'Homo sapiens': '', 'reads': 197}}}, 'reads': 197}}, 'reads': 197}, 'reads': 197}

    @patch('analysis_driver.quality_control.ContaminationBlast.get_rank')
    @patch('analysis_driver.quality_control.ContaminationBlast.ncbi')
    @patch('analysis_driver.quality_control.ContaminationBlast.taxid_translator')
    def test_get_all_taxa_identified2(self, mocked_taxid_translator, mocked_ncbi, mocked_rank):
        mocked_rank.return_value = {99802: 'species',
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
                                    1224679: 'order'}

        mocked_ncbi.return_value = None
        mocked_taxid_translator.return_value = ('Metazoa', 'Platyhelminthes', 'Cestoda', 'Diphyllobothriidea', 'Diphyllobothriidae', 'Spirometra', 'Spirometra erinaceieuropaei')
        taxon_dict = {'Metazoa': {'Chordata': {'reads': 197, 'Mammalia': {'Primates': {'reads': 197, 'Hominidae': {'reads': 197, 'Homo': {'Homo sapiens': '', 'reads': 197}}}, 'reads': 197}}, 'reads': 197}, 'reads': 197}
        taxon = '99802'
        taxids = Counter({'9606': 197, '99802': 2, '9598': 2})
        test_spirometra_and_human = self.contamination_blast.get_all_taxa_identified(taxon_dict, taxon, taxids)
        assert test_spirometra_and_human == {'reads': 199, 'Metazoa': {'Chordata': {'reads': 197, 'Mammalia': {'reads': 197, 'Primates': {'Hominidae': {'Homo': {'reads': 197, 'Homo sapiens': ''}, 'reads': 197}, 'reads': 197}}}, 'reads': 2, 'Platyhelminthes': {'Cestoda': {'Diphyllobothriidea': {'reads': 2, 'Diphyllobothriidae': {'reads': 2, 'Spirometra': {'reads': 2, 'Spirometra erinaceieuropaei': ''}}}, 'reads': 2}, 'reads': 2}}}

    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_sample_fastq(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        sample_fastq_outfile = self.contamination_blast.run_sample_fastq(self.fastq_file)
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
        blast_outfile = self.contamination_blast.run_blast(self.fasta_file)
        assert blast_outfile == 'test_run/fastaFile1_blastn'
        mocked_execute.assert_called_once_with(
            "export PATH=$PATH:/path/to/db/dir; path/to/blastn -query fastaFile1.fasta -db path/to/db/dir/nt -out test_run/fastaFile1_blastn -num_threads 12 -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'",
            working_dir='test_run',
            mem=20,
            cpus=12,
            job_name='contamination_blast'
        )
