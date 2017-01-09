import os.path
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
        self.fastq_file = fastq_file
        self._ncbi = None

    def sample_fastq_command(self, fastq_file, nb_reads):
        fasta_outfile = os.path.join(
            self.working_dir,
            os.path.basename(fastq_file).split('.')[0] + '_sample%s.fasta' % nb_reads
        )
        cmd = 'set -o pipefail; {seqtk} sample {fastq} {nb_reads} | {seqtk} seq -a > {fasta}'.format(
            nb_reads=nb_reads, seqtk=cfg['tools']['seqtk'], fastq=fastq_file, fasta=fasta_outfile
        )
        return cmd, fasta_outfile

    def fasta_blast_command(self, fasta_file):
        db_dir = cfg['contamination-check']['db_dir']
        nt_db = os.path.join(db_dir, 'nt')
        blast_outfile = os.path.join(self.working_dir, os.path.basename(fasta_file).split('.')[0] + '_blastn')

        cmd = ('export PATH=$PATH:/{db_dir}; {blastn} -query {fasta_file} -db {nt} -out {blast_outfile} '
               '-num_threads 12 -max_target_seqs 1 -max_hsps 1 '
               "-outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'")
        cmd = cmd.format(db_dir=db_dir, blastn=cfg['tools']['blastn'], fasta_file=fasta_file, nt=nt_db,
                         blast_outfile=blast_outfile)
        return cmd, blast_outfile

    @staticmethod
    def get_taxids(blast):
        taxids = Counter()
        with open(blast) as f:
            for line in f:
                taxid = line.split()[7]
                # sometime more than one taxid are reported for a specific hit
                # they're all resolving to the same tax name
                taxid = taxid.split(';')[0]
                taxids[taxid] += 1
        return taxids

    @property
    def ncbi(self):
        if not self._ncbi:
            db_path = cfg['contamination-check']['ete_db']
            if os.path.exists(db_path):
                self._ncbi = NCBITaxa(dbfile=db_path)
            else:
                raise AnalysisDriverError('Cannot locate the ETE taxon database')
        return self._ncbi

    def get_ranks(self, taxon):
        """Retrieve the rank of each of the taxa from that taxid's lineage"""
        if taxon != 'N/A':
            try:
                l = self.ncbi.get_lineage(int(taxon))
                rank = self.ncbi.get_rank(l)
            except ValueError:
                rank = {0: 'rank unavailable'}
                self.warning('The taxid %s does not exist in the ETE TAXDB' % taxon)
            return rank

    def get_all_taxa_identified(self, taxon_dict, taxon, taxids):
        num_reads = taxids[taxon]
        ranks = self.get_ranks(taxon)
        required_ranks = ['superkingdom', 'kingdom', 'phylum', 'class',  'order', 'family', 'genus', 'species']

        taxon_dict_for_current_rank = taxon_dict
        for required_rank in required_ranks:
            taxid_for_rank = [i for i in ranks if ranks[i] == required_rank]
            # list of one taxid for that rank because get_taxid_translator requires a list
            if taxid_for_rank:
                taxon_for_rank = list(self.ncbi.get_taxid_translator(taxid_for_rank).values()).pop()
            else:
                taxon_for_rank = 'Unavailable'
            if taxon_for_rank and taxon_for_rank not in taxon_dict_for_current_rank:
                taxon_dict_for_current_rank[taxon_for_rank] = {'reads': num_reads}
            elif taxon_for_rank:
                taxon_dict_for_current_rank[taxon_for_rank]['reads'] += num_reads
            taxon_dict_for_current_rank = taxon_dict_for_current_rank[taxon_for_rank]

        return taxon_dict

    def run_sample_fastq(self, fastq_file, nb_reads):
        self.dataset.start_stage('sample_fastq')
        sample_fastq_command, fasta_outfile = self.sample_fastq_command(fastq_file, nb_reads)
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
        nb_reads = 3000
        fasta_outfile = self.run_sample_fastq(self.fastq_file[0], nb_reads)
        blast_outfile = self.run_blast(fasta_outfile)
        taxids = self.get_taxids(blast_outfile)
        taxon_dict = {'Total': nb_reads}
        for taxon in taxids:
            taxon_dict = self.get_all_taxa_identified(taxon_dict, taxon, taxids)
        outpath = os.path.join(self.working_dir, 'taxa_identified.json')
        with open(outpath, 'w') as outfile:
            json.dump(taxon_dict, outfile, sort_keys=True, indent=4, separators=(',', ':'))

    def run(self):
        try:
            self.check_for_contamination()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
