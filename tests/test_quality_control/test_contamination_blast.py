import os
from tests.test_quality_control.qc_tester import QCTester
from analysis_driver.quality_control import ContaminationBlast
from unittest.mock import patch

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
        assert command == 'path/to/seqtk sample fastqFile1.fastq 3000 | seqtk seq -a > test_run/test_run/fastqFile1_sample3000.fasta'
        assert outfile == 'test_run/test_run/fastqFile1_sample3000.fasta'

    def test_fasta_blast_command(self):
        command, outfile = self.contamination_blast.fasta_blast_command(self.fasta_file)
        assert command == "path/to/blastn -query fastaFile1.fasta -db path/to/nt_db -out fastaFile1_blastn -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"
        assert outfile =='fastaFile1_blastn'

    def test_retrieve_identified_taxa(self):
        blast_file = os.path.join(self.assets_path, 'blast_outfile')
        taxa_retrieved = self.contamination_blast.retrieve_identified_taxa(blast_file)
        assert taxa_retrieved == {'Hominidae': 199, 'Diphyllobothriidae': 2}

    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_sample_fastq(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        sample_fastq_outfile = self.contamination_blast.run_sample_fastq(self.fastq_file)
        assert sample_fastq_outfile == 'test_run/test_run/fastqFile1_sample3000.fasta'
        mocked_execute.assert_called_once_with(
            'path/to/seqtk sample fastqFile1.fastq 3000 | seqtk seq -a > test_run/test_run/fastqFile1_sample3000.fasta',
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
        assert blast_outfile == 'fastaFile1_blastn'
        mocked_execute.assert_called_once_with(
            "path/to/blastn -query fastaFile1.fasta -db path/to/nt_db -out fastaFile1_blastn -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'",
            working_dir='test_run',
            mem=10,
            cpus=2,
            job_name='contamination_blast'
        )
