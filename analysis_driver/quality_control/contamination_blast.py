import os
import json
from ete3 import NCBITaxa
from collections import Counter
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
        db_dir = cfg['contamination-check']['db_dir']
        nt_db = os.path.join(db_dir, 'nt')
        blast_outfile = fasta_file.split('.')[0] + '_blastn'
        blast_cmd = "export PATH=$PATH:/%s; %s -query %s -db %s -out %s -num_threads 12 -outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'" % (db_dir, blastn_bin, fasta_file, nt_db, blast_outfile)
        return blast_cmd, blast_outfile

    def check_db_exists(self, db_path):
        if os.path.exists(db_path):
            return True
        else:
            raise AnalysisDriverError('Cannot locate the ETE taxon database')

    def get_taxids(self, blast):
        taxids = Counter()
        with open(blast) as openfile:
            blast = openfile.readlines()
            for line in blast:
                taxid = line.split()[7]
                taxids[taxid] +=1
        return taxids

    def get_all_taxa_identified(self, taxon_dict, taxon, taxids, ncbi):
        num_reads = taxids[taxon]
        l = ncbi.get_lineage(taxon)
        rank = ncbi.get_rank(l)
        kingdom_taxid = [i for i in rank if rank[i] == 'kingdom']
        phylum_taxid = [i for i in rank if rank[i] == 'phylum']
        class_taxid = [i for i in rank if rank[i] == 'class']
        order_taxid = [i for i in rank if rank[i] == 'order']
        family_taxid = [i for i in rank if rank[i] == 'family']
        genus_taxid = [i for i in rank if rank[i] == 'genus']
        species_taxid = [i for i in rank if rank[i] == 'species']
        taxon_kingdom = list((ncbi.get_taxid_translator(kingdom_taxid)).values())[0]
        taxon_phylum = list((ncbi.get_taxid_translator(phylum_taxid)).values())[0]
        taxon_class = list((ncbi.get_taxid_translator(class_taxid)).values())[0]
        taxon_order = list((ncbi.get_taxid_translator(order_taxid)).values())[0]
        taxon_family = list((ncbi.get_taxid_translator(family_taxid)).values())[0]
        taxon_genus = list((ncbi.get_taxid_translator(genus_taxid)).values())[0]
        taxon_species = list((ncbi.get_taxid_translator(species_taxid)).values())[0]
        if taxon_kingdom not in taxon_dict:
            taxon_dict[taxon_kingdom] = {}
            taxon_dict['reads'] = num_reads
        else:
            taxon_dict['reads'] += num_reads
        if taxon_phylum not in taxon_dict.get(taxon_kingdom, {}):
            taxon_dict[taxon_kingdom][taxon_phylum] = {}
            taxon_dict[taxon_kingdom]['reads'] = num_reads
        else:
            taxon_dict[taxon_kingdom]['reads'] += num_reads
        if taxon_class not in taxon_dict.get(taxon_kingdom, {}).get(taxon_phylum, {}):
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class] = {}
            taxon_dict[taxon_kingdom][taxon_phylum]['reads'] = num_reads
        else:
            taxon_dict[taxon_kingdom][taxon_phylum]['reads'] += num_reads
        if taxon_order not in taxon_dict.get(taxon_kingdom, {}).get(taxon_phylum, {}).get(taxon_class, {}):
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order] = {}
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class]['reads'] = num_reads
        else:
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class]['reads'] += num_reads
        if taxon_family not in taxon_dict.get(taxon_kingdom, {}).get(taxon_phylum, {}).get(taxon_class, {}).get(taxon_order, {}):
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family] = {}
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order]['reads'] = num_reads
        else:
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order]['reads'] += num_reads
        if taxon_genus not in taxon_dict.get(taxon_kingdom, {}).get(taxon_phylum, {}).get(taxon_class, {}).get(taxon_order, {}).get(taxon_family, {}):
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family][taxon_genus] = {}
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family]['reads'] = num_reads
        else:
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family]['reads'] += num_reads
        if taxon_species not in taxon_dict.get(taxon_kingdom, {}).get(taxon_phylum, {}).get(taxon_class, {}).get(taxon_order, {}).get(taxon_family, {}).get(taxon_genus, {}):
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family][taxon_genus][taxon_species] = ''
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family][taxon_genus]['reads'] = num_reads
        else:
            taxon_dict[taxon_kingdom][taxon_phylum][taxon_class][taxon_order][taxon_family][taxon_genus]['reads'] += num_reads

        return taxon_dict

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
            mem=20
        )
        exit_status = contamination_blast_executor.join()
        self.dataset.end_stage('contamination_blast', exit_status)
        return blast_outfile

    def check_for_contamination(self):
        fasta_outfile = self.run_sample_fastq(self.fastq_file)
        blast_outfile = self.run_blast(fasta_outfile)
        taxids = self.get_taxids(blast_outfile)
        db_path = cfg['contamination-check']['ete_db']
        if self.check_db_exists(db_path):
            ncbi = NCBITaxa(dbfile=db_path)
            taxon_dict = {}
            for taxon in taxids:
                taxon_dict = self.get_all_taxa_identified(taxon_dict, taxon, taxids, ncbi)


        outpath = os.path.join(self.working_dir, 'taxa_identified.json')
        with open(outpath, 'w') as outfile:
            taxa_identified_json = json.dumps(taxon_dict,
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