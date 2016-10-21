import os
import json
from ete3 import NCBITaxa
from egcg_core import executor
from .quality_control_base import QualityControl
from analysis_driver.config import default as cfg
from analysis_driver.exceptions import AnalysisDriverError

class ContaminationBlast(QualityControl):

    def __init__(self, dataset, working_dir, fastq_file):
        super().__init__(dataset, working_dir)
        self.working_dir = working_dir
        self.fastq_file = fastq_file


    def sample_fastq_command(self, fastq_file):
        seqtk_bin = cfg['tools']['seqtk']
        fasta_outfile = os.path.join(self.working_dir, self.dataset.name, fastq_file.split('.')[0] + '_sample3000.fasta')
        seqtk_sample_cmd = 'set -o pipefail; {seqtk} sample {fastq} 3000 | {seqtk} seq -a > {fasta}'.format(seqtk=seqtk_bin, fastq=fastq_file, fasta=fasta_outfile)
        return seqtk_sample_cmd, fasta_outfile


    def fasta_blast_command(self, fasta_file):
        blastn_bin = cfg['tools']['blastn']
        nt_db = cfg['contamination-check']['nt_db']
        blast_outfile = fasta_file.split('.')[0] + '_blastn'
        blast_cmd = "%s -query %s -db %s -out %s -num_threads 12 -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'" % (blastn_bin, fasta_file, nt_db, blast_outfile)
        return blast_cmd, blast_outfile

    def check_db_exists(self, db_path):
        if os.path.exists(db_path):
            return True
        else:
            raise AnalysisDriverError('Cannot locate the ETE taxon database')


    def taxon_info(self, taxon, db_path):
        if self.check_db_exists(db_path):
            ncbi = NCBITaxa(dbfile=db_path)
            lineage = ncbi.get_lineage(taxon)
            names = ncbi.get_taxid_translator(lineage)
            rank = ncbi.get_rank(lineage)
            class_taxid = [i for i in rank if rank[i] == 'family']
            return names, class_taxid


    def get_taxids(self, blast):
        taxids = {}
        with open(blast) as openfile:
            blast = openfile.readlines()
            for line in blast:
                taxid = line.split()[7]
                if taxid not in taxids:
                    taxids[taxid] = 1
                else:
                    taxids[taxid] += 1
        return taxids


    def get_all_taxa_identified(self, taxon, taxids, db_path):
        taxa_identified = {}
        num_reads = taxids.get(taxon)
        names, class_taxid = self.taxon_info(taxon, db_path)
        if class_taxid:
            taxon_id = names.get(class_taxid[0])
            if taxon_id not in taxa_identified:
                taxa_identified[taxon_id] = num_reads
            else:
                taxa_identified[taxon_id]+=num_reads
        return taxa_identified


    def retrieve_identified_taxa(self, blast):
        taxa_identified = {}
        db_path = cfg['contamination-check']['ete_db']
        taxids = self.get_taxids(blast)
        for taxon in taxids:
            taxa_identified = self.get_all_taxa_identified(taxon, taxids, db_path)
        return taxa_identified


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
            cpus=12,
            mem=10
        )
        exit_status = contamination_blast_executor.join()
        self.dataset.end_stage('contamination_blast', exit_status)
        return blast_outfile


    def check_for_contamination(self):
        taxa_identified = {}
        fasta_outfile = self.run_sample_fastq(self.fastq_file)
        blast_outfile = self.run_blast(fasta_outfile)
        taxa = self.retrieve_identified_taxa(blast_outfile)
        for taxon in taxa:
            if taxon not in taxa_identified:
                taxa_identified[taxon] = taxa[taxon]
            else:
                taxa_identified[taxon] += taxa[taxon]

        outpath = os.path.join(self.working_dir, 'taxa_identified.json')
        with open(outpath, 'w') as outfile:
            taxa_identified_json = json.dumps(taxa_identified,
                                            sort_keys=True, indent=4,
                                            separators=(',', ':'))
            outfile.write(taxa_identified_json)


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