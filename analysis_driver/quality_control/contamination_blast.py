import os
from ete3 import NCBITaxa
from egcg_core import executor
from .quality_control_base import QualityControl
from analysis_driver.config import default as cfg

class ContaminationBlast(QualityControl):
    def __init__(self, dataset, working_dir, fastq_files):
        super().__init__(dataset, working_dir)
        self.working_dir = working_dir
        self.fastq_files = fastq_files

    def sample_fastq_command(self, fastq_file):
        seqtk_bin = cfg['tools']['seqtk']
        fasta_outfile = os.path.join(self.working_dir, self.dataset.name, fastq_file.split('.')[0] + '_sample3000.fasta')
        seqtk_sample_cmd = '%s sample %s 3000 | seqtk seq -a > %s' % (seqtk_bin, fastq_file, fasta_outfile)
        return seqtk_sample_cmd, fasta_outfile

    def fasta_blast_command(self, fasta_file):
        blastn_bin = cfg['tools']['blastn']
        nt_db = cfg['contamination-check']['nt_db']
        blast_outfile = fasta_file.split('.')[0] + '_blastn'
        blast_cmd = "%s -query %s -db %s -out %s -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'" % (blastn_bin, fasta_file, nt_db, blast_outfile)
        return blast_cmd, blast_outfile

    def retrieve_identified_taxa(self, blast):
        taxids = {}
        classes = {}
        ncbi = NCBITaxa()

        with open(blast) as openfile:
            blast = openfile.readlines()
            for line in blast:
                taxid = line.split()[7]
                if taxid not in taxids:
                    taxids[taxid] = 1
                else:
                    taxids[taxid] += 1

        for taxon in taxids:
            num_reads = taxids[taxon]
            lineage = ncbi.get_lineage(taxon)
            names = ncbi.get_taxid_translator(lineage)
            rank = ncbi.get_rank(lineage)
            class_taxid = [i for i in rank if rank[i] == 'family']
            if class_taxid:
                class_name = names[class_taxid[0]]
                if class_name not in classes:
                    classes[class_name] = num_reads
                else:
                    classes[class_name]+=num_reads
        return classes

    def run_sample_fastq(self, fastq_file):
        self.dataset.start_stage('sample_fastq')
        sample_fastq_command, fasta_outfile = self.sample_fastq_command(fastq_file)
        sample_fastq_executor = executor.execute(
            sample_fastq_command,
            job_name='sample_fastq',
            working_dir=self.working_dir,
            cpus=2,
            mem=10
        )
        exit_status = sample_fastq_executor.join()
        self.dataset.end_stage('sample_fastq', exit_status)
        return fasta_outfile

    def run_blast(self, fasta_file):
        self.dataset.start_stage('contamination_blast')
        fasta_blast_command, blast_outfile = self.fasta_blast_command(fasta_file)
        contamination_blast_executor = executor.execute(
            fasta_blast_command,
            job_name='contamination_blast',
            working_dir=self.working_dir,
            cpus=2,
            mem=10
        )
        exit_status = contamination_blast_executor.join()
        self.dataset.end_stage('contamination_blast', exit_status)
        return blast_outfile

    def check_for_contamination(self):
        taxa_identified = {}
        for fastq_file in self.fastq_files:

            fasta_outfile = self.run_sample_fastq(fastq_file)
            blast_outfile = self.run_blast(fasta_outfile)
            taxa = self.retrieve_identified_taxa(blast_outfile)

            for taxon in taxa:
                if taxon not in taxa_identified:
                    taxa_identified[taxon] = taxa[taxon]
                else:
                    taxa_identified[taxon] += taxa[taxon]
        return taxa_identified

    def run(self):
        try:
            self.taxa_identified = self.check_for_contamination()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.taxa_identified