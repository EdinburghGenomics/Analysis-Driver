import os
import shutil
import json
import os.path
from collections import Counter
from ete3 import NCBITaxa
from luigi import Parameter, IntParameter
from egcg_core import executor, util
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg
from analysis_driver.segmentation import Stage
from analysis_driver.tool_versioning import toolset


class FastqScreen(Stage):
    fq_pattern = Parameter()

    def _fastqscreen_command(self):
        fqs = util.find_files(self.fq_pattern)
        assert 1 <= len(fqs) <= 2, 'Bad number of fastqs: %s' % fqs
        return '%s --aligner bowtie2 %s --conf %s --force' % (
            toolset['fastqscreen'], ' '.join(fqs), cfg['contamination-check']['fastqscreen_conf']
        )

    @property
    def fastqscreen_expected_outfiles(self):
        return self.fq_pattern.replace('.fastq.gz', '_screen.txt')

    def _run(self):
        fastqscreen_run_command = self._fastqscreen_command()
        return executor.execute(
            fastqscreen_run_command,
            job_name='fastqscreen',
            working_dir=self.job_dir,
            cpus=2,
            mem=10
        ).join()


class VerifyBamID(Stage):
    bam_file = Parameter()

    @property
    def filtered_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '_chr22.bam')

    def _filter_bam(self):
        # use only chromosome 22 for speed
        return executor.execute(
            toolset['samtools'] + ' view -b %s chr22 > %s' % (util.find_file(self.bam_file), self.filtered_bam),
            job_name='filter_bam22',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()

    def _index_filtered_bam(self):
        return executor.execute(
            toolset['samtools'] + ' index %s' % self.filtered_bam,
            job_name='index_bam22',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()

    def _verify_bam_id(self):
        cmd = '%s --bam %s --vcf %s --out %s' % (
            toolset['verifybamid'],
            self.filtered_bam,
            cfg['contamination-check']['population_vcf'],
            os.path.join(self.job_dir, self.dataset.name + '-chr22-vbi')
        )
        exit_status = executor.execute(
            cmd,
            job_name='verify_bam_id',
            working_dir=self.job_dir,
            cpus=1,
            mem=4
        ).join()
        sample_vbi_self = os.path.join(self.job_dir, self.dataset.name + '-chr22-vbi.selfSM')
        if os.path.exists(sample_vbi_self):
            bam_dir = os.path.dirname(self.bam_file)
            dest = os.path.join(bam_dir, os.path.basename(sample_vbi_self))
            shutil.copyfile(sample_vbi_self, dest)
        return exit_status

    def _run(self):
        exit_status = 0
        exit_status += self._filter_bam()
        exit_status += self._index_filtered_bam()
        exit_status += self._verify_bam_id()
        return exit_status


class VCFStats(Stage):
    vcf_file = Parameter()

    def _run(self):
        vcf = util.find_file(self.vcf_file)
        if not vcf:
            return 1

        name, ext = os.path.splitext(vcf)
        stats_file = name + '.stats'
        cmd = '%s vcfstats %s > %s' % (toolset['rtg'], vcf, stats_file)
        exit_status = executor.execute(
            cmd,
            job_name='rtg_vcfstats',
            working_dir=self.job_dir,
            cpus=4,
            mem=32
        ).join()
        return exit_status


class Blast(Stage):
    fastq_file = Parameter()
    nb_reads = IntParameter(default=3000)
    _ncbi = None

    @property
    def fasta_outfile(self):
        return os.path.join(
            self.job_dir,
            os.path.basename(self.fastq_file).split('.')[0] + '_sample%s.fasta' % self.nb_reads
        )

    @property
    def blast_outfile(self):
        return os.path.join(self.job_dir, os.path.basename(self.fasta_outfile).split('.')[0] + '_blastn')

    def sample_fastq_command(self):
        return 'set -o pipefail; {seqtk} sample {fastq} {nb_reads} | {seqtk} seq -a > {fasta}'.format(
            nb_reads=self.nb_reads, seqtk=toolset['seqtk'], fastq=util.find_file(self.fastq_file),
            fasta=self.fasta_outfile
        )

    def fasta_blast_command(self):
        db_dir = cfg['contamination-check']['db_dir']
        cmd = ('export PATH=$PATH:/{db_dir}; {blastn} -query {fasta_file} -db {nt} -out {blast_outfile} '
               '-num_threads 12 -max_target_seqs 1 -max_hsps 1 '
               "-outfmt '6 qseqid sseqid length pident evalue sgi sacc staxids sscinames scomnames stitle'")
        return cmd.format(
            db_dir=db_dir, blastn=toolset['blastn'], fasta_file=self.fasta_outfile,
            nt=os.path.join(db_dir, 'nt'), blast_outfile=self.blast_outfile
        )

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

    def run_sample_fastq(self):
        return executor.execute(
            self.sample_fastq_command(),
            job_name='sample_fastq',
            working_dir=self.job_dir,
            cpus=2,
            mem=10
        ).join()

    def run_blast(self):
        return executor.execute(
            self.fasta_blast_command(),
            job_name='contamination_blast',
            working_dir=self.job_dir,
            cpus=12,
            mem=20
        ).join()

    def _run(self):
        exit_status = self.run_sample_fastq()
        exit_status += self.run_blast()
        taxids = self.get_taxids(self.blast_outfile)
        taxon_dict = {'Total': self.nb_reads}
        for taxon in taxids:
            taxon_dict = self.get_all_taxa_identified(taxon_dict, taxon, taxids)

        outpath = os.path.join(self.job_dir, 'taxa_identified.json')
        with open(outpath, 'w') as outfile:
            json.dump(taxon_dict, outfile, sort_keys=True, indent=4, separators=(',', ':'))

        return exit_status
