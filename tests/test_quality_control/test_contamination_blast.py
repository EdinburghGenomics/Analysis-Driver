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
        assert command == 'set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | path/to/seqtk seq -a > test_run/test_run/fastqFile1_sample3000.fasta'
        assert outfile == 'test_run/test_run/fastqFile1_sample3000.fasta'

    def test_fasta_blast_command(self):
        command, outfile = self.contamination_blast.fasta_blast_command(self.fasta_file)
        assert command == "path/to/blastn -query fastaFile1.fasta -db path/to/nt_db -out fastaFile1_blastn -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'"
        assert outfile =='fastaFile1_blastn'

    def test_get_taxids(self):
        blast_file = os.path.join(self.assets_path, 'blast_outfile')
        taxids = self.contamination_blast.get_taxids(blast_file)
        assert taxids == {'9598': 2, '9606': 197, '99802': 2}

    def test_update_classes(self):
        taxon = '9606'
        classes = {}
        taxids = {'99802': 2, '9598': 2, '9606': 197}
        db_path = 'test/db/path'
        patch_names = {1: u'root',
                       33154: u'Opisthokonta',
                       9347: u'Eutheria',
                       9604: u'Hominidae',
                       9605: u'Homo',
                       9606: u'Homo sapiens',
                       1338369: u'Dipnotetrapodomorpha',
                       32523: u'Tetrapoda',
                       32524: u'Amniota',
                       32525: u'Theria',
                       376913: u'Haplorrhini',
                       7711: u'Chordata',
                       314146: u'Euarchontoglires',
                       314293: u'Simiiformes',
                       9526: u'Catarrhini',
                       314295: u'Hominoidea',
                       33208: u'Metazoa',
                       33213: u'Bilateria',
                       7742: u'Vertebrata',
                       117570: u'Teleostomi',
                       117571: u'Euteleostomi',
                       2759: u'Eukaryota',
                       6072: u'Eumetazoa',
                       1437010: u'Boreoeutheria',
                       8287: u'Sarcopterygii',
                       7776: u'Gnathostomata',
                       40674: u'Mammalia',
                       9443: u'Primates',
                       33511: u'Deuterostomia',
                       207598: u'Homininae',
                       131567: u'cellular organisms',
                       89593: u'Craniata'}
        patch_class_taxid = [9604]
        with patch('analysis_driver.quality_control.ContaminationBlast.taxon_info', return_value=(patch_names, patch_class_taxid)):
            classes = self.contamination_blast.update_classes(taxon, classes, taxids, db_path)
            assert classes == {'Hominidae': 197}

        classes = {'Hominidae': 197}
        with patch('analysis_driver.quality_control.ContaminationBlast.taxon_info', return_value=(patch_names, patch_class_taxid)):
            classes = self.contamination_blast.update_classes(taxon, classes, taxids, db_path)
            assert classes == {'Hominidae': 394}

        classes = {'Diphyllobothriidae': 3}
        with patch('analysis_driver.quality_control.ContaminationBlast.taxon_info', return_value=(patch_names, patch_class_taxid)):
            classes = self.contamination_blast.update_classes(taxon, classes, taxids, db_path)
            assert classes == {'Diphyllobothriidae': 3, 'Hominidae': 197}


    @patch('analysis_driver.dataset.rest_communication')
    @patch('egcg_core.executor.execute')
    def test_run_sample_fastq(self, mocked_execute, mocked_rest):
        instance = mocked_execute.return_value
        instance.join.return_value = 0
        sample_fastq_outfile = self.contamination_blast.run_sample_fastq(self.fastq_file)
        assert sample_fastq_outfile == 'test_run/test_run/fastqFile1_sample3000.fasta'
        mocked_execute.assert_called_once_with(
            'set -o pipefail; path/to/seqtk sample fastqFile1.fastq 3000 | path/to/seqtk seq -a > test_run/test_run/fastqFile1_sample3000.fasta',
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
            cpus=12,
            job_name='contamination_blast'
        )
