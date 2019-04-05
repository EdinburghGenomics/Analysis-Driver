import os
import shutil
from os.path import join
from unittest.mock import patch, PropertyMock, MagicMock

import luigi
import sys
from egcg_core import executor, util
from egcg_core.config import cfg
from egcg_core.constants import ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_ELEMENT_ID, ELEMENT_PROJECT_ID, \
    ELEMENT_RUN_NAME, ELEMENT_LANE

from analysis_driver import segmentation
from analysis_driver.config import load_config
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver import quality_control as qc
from analysis_driver.pipelines.common import bgzip_and_tabix, tabix_vcf, MergeFastqs, SamtoolsStats
from analysis_driver.segmentation import Parameter
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.bash_commands import picard_command

toolset_type = 'gatk4_sample_processing'
name = 'variant_calling_gatk4'


class SplitFastqStage(segmentation.Stage):
    """
    Generic stage that provides functionality to access
       - fastq files for each run element
       - fastq files for each run element recompressed and indexed
       - chunks used to access the indexed fastq files
    """

    def _find_fastqs_for_run_element(self, run_element):
        local_fastq_dir = os.path.join(cfg['sample']['input_dir'], run_element.get(ELEMENT_RUN_NAME))
        self.debug('Searching for fastqs in ' + local_fastq_dir)
        return util.find_fastqs(
            local_fastq_dir,
            run_element.get(ELEMENT_PROJECT_ID),
            self.dataset.name,
            run_element.get(ELEMENT_LANE)
        )

    @property
    def indexed_fastq_dir(self):
        return os.path.join(self.job_dir, 'fastq_indexed')

    def _indexed_fastq_for_run_element(self, run_element):
        return [
            os.path.join(self.indexed_fastq_dir, run_element.get(ELEMENT_RUN_ELEMENT_ID) + fastq_file[-len('_R1_001.fastq.gz'):])
            for fastq_file in self._find_fastqs_for_run_element(run_element)
        ]

    @staticmethod
    def chunks_from_fastq(indexed_fq_files):
        split_lines = 100000000
        grabix_index = indexed_fq_files[0] + '.gbi'
        assert os.path.isfile(grabix_index)
        with open(grabix_index) as open_file:
            next(open_file)  # discard first line of the file
            nb_lines = int(next(open_file))  # second line contains the number of lines in the file
        chunks = []
        last = 1
        for chunki in range(nb_lines // split_lines + min(1, nb_lines % split_lines)):
            new = last + split_lines - 1
            chunks.append((last, min(new, nb_lines)))
            last = new + 1
        return chunks

    @property
    def split_alignment_dir(self):
        return join(self.job_dir, 'split_alignment')

    def chuncked_bam_file(self, run_element, chunk):
        return join(self.split_alignment_dir,
                    run_element.get(ELEMENT_RUN_ELEMENT_ID) + '_name_sort_%s_%s' % chunk + '.bam')


class FastqIndex(SplitFastqStage):
    """get all fastq files, uncompress and recompress, then index with grabix."""

    def _run(self):
        os.makedirs(self.indexed_fastq_dir, exist_ok=True)
        os.makedirs(self.split_alignment_dir, exist_ok=True)
        commands = []
        for run_element in self.dataset.run_elements:
            if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                for ifastqs, ofastqs in zip(
                        self._find_fastqs_for_run_element(run_element),
                        self._indexed_fastq_for_run_element(run_element)
                ):
                    cmd_template = 'gunzip -c {ipath} | {pbgzip} -n 16  -c /dev/stdin > {opath}'
                    commands.append(cmd_template.format(ipath=ifastqs, pbgzip=toolset['pbgzip'], opath=ofastqs))
        exit_status = executor.execute(
            *commands,
            job_name='compress_fastq',
            working_dir=self.job_dir,
            log_commands=False
        ).join()

        if exit_status == 0:
            commands = []
            for run_element in self.dataset.run_elements:
                if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                    for fastq_file in self._indexed_fastq_for_run_element(run_element):
                        commands.append(toolset['grabix'] + ' index ' + fastq_file)
            exit_status = executor.execute(
                *commands,
                job_name='index_fastq',
                working_dir=self.job_dir
            ).join()
        return exit_status


class SplitBWA(SplitFastqStage):
    """
    Run bwa on check of fastq file provided by SplitFastqStage.
    """

    def bwa_command(self, fastq_pair, reference, expected_output_bam, read_group, chunk):
        tmp_file = expected_output_bam
        read_group_str = r'@RG\t%s' % r'\t'.join(['%s:%s' % (k, read_group[k]) for k in sorted(read_group)])

        command_bwa = '{bwa} mem -K 100000000 -Y -R \'{read_group}\' -M -t 2 {ref} ' \
                      '<({grabix} grab {fastq1} {chunk}) ' \
                      '<({grabix} grab {fastq2} {chunk})'
        command_bwa = command_bwa.format(bwa=toolset['bwa'], read_group=read_group_str, ref=reference,
                                         grabix=toolset['grabix'], fastq1=fastq_pair[0], fastq2=fastq_pair[1],
                                         chunk='%s %s' % (chunk[0], chunk[1]))

        command_samtools = '{samtools} sort -n -m 1G -O bam -T {tmp_file} -o {bam_file} -'
        command_samtools = command_samtools.format(samtools=toolset['samtools'], tmp_file=tmp_file,
                                                   bam_file=expected_output_bam)

        cmd = 'set -o pipefail; ' + ' | '.join([command_bwa, command_samtools])
        return cmd

    def _run(self):
        """Get all fastq files, and align them in chunks."""

        commands = []
        for run_element in self.dataset.run_elements:
            if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                indexed_fq_files = self._indexed_fastq_for_run_element(run_element)
                for chunk in self.chunks_from_fastq(indexed_fq_files):
                    commands.append(self.bwa_command(
                        fastq_pair=indexed_fq_files,
                        reference=self.dataset.reference_genome,
                        expected_output_bam=self.chuncked_bam_file(run_element, chunk),
                        read_group={'ID': run_element[ELEMENT_RUN_ELEMENT_ID], 'PU': run_element[ELEMENT_RUN_ELEMENT_ID],
                                    'SM': self.dataset.user_sample_id, 'PL': 'illumina'},
                        chunk=chunk
                    ))
        exit_status = executor.execute(
            *commands,
            job_name='bwa_split',
            working_dir=self.job_dir,
            cpus=2,
            mem=12
        ).join()
        return exit_status


class MergeBamAndDup(SplitFastqStage):
    """
    Merge bam file generated in bwa mem and defined in SplitFastqStage.
    """

    def all_chunk_bam_list_file(self):
        all_chunk_bams = []
        for run_element in self.dataset.run_elements:
            indexed_fq_files = self._indexed_fastq_for_run_element(run_element)
            for chunk in self.chunks_from_fastq(indexed_fq_files):
                all_chunk_bams.append(self.chuncked_bam_file(run_element, chunk))
        bam_list_file = os.path.join(self.job_dir, self.dataset.name + '_all_bam_files.list')
        with open(bam_list_file, 'w') as open_file:
            open_file.write('\n'.join(all_chunk_bams))
        return bam_list_file

    @property
    def cat_tmp_dir(self):
        return os.path.join(self.job_dir, 'cat_tmp')

    @property
    def dup_tmp_dir(self):
        return os.path.join(self.job_dir, 'dup_tmp')

    @property
    def sorted_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')

    def merge_command(self):
        cat_cmd = '{bamcat} level=0 tmpfile={cat_tmp} `cat {bam_list_file}`'.format(
            bamcat=toolset['bamcat'],
            cat_tmp=os.path.join(self.cat_tmp_dir, self.dataset.name),
            bam_list_file=self.all_chunk_bam_list_file()
        )
        bamsormadup_cmd = '{bamsormadup} threads=16 tmpfile={dup_tmp} indexfilename={merged_bam}.bai > {merged_bam}'.format(
            bamsormadup=toolset['biobambam_sortmapdup'],
            dup_tmp=os.path.join(self.dup_tmp_dir, self.dataset.name),
            merged_bam=self.sorted_bam
        )
        return 'set -o pipefail; ' + ' | '.join([cat_cmd, bamsormadup_cmd])

    def _run(self):
        os.makedirs(self.cat_tmp_dir)
        os.makedirs(self.dup_tmp_dir)
        exit_status = executor.execute(
            self.merge_command(),
            job_name='merge_dup_bam',
            working_dir=self.job_dir,
            cpus=6,
            mem=36
        ).join()
        shutil.rmtree(self.cat_tmp_dir)
        shutil.rmtree(self.dup_tmp_dir)
        return exit_status


class GATK4Stage(segmentation.Stage):

    @property
    def exec_dir(self):
        d = os.path.join(self.job_dir, 'slurm_and_logs')
        os.makedirs(d, exist_ok=True)
        return d

    @property
    def gatk_run_dir(self):
        d = os.path.join(self.job_dir, 'gatk4')
        os.makedirs(d, exist_ok=True)
        return d

    def _gatk_cmd(self, run_cls, output=None, memory=2, input=None, spark_core=1, reference='default', ext=None,
                  picard_tool=False):
        base_cmd = toolset['gatk4_bin'] + \
                   ' --java-options "-Djava.io.tmpdir={tmp_dir} -XX:+UseSerialGC -Xmx{memory}G" ' + \
                   '{run_cls} '
        if picard_tool:
            input_param = 'INPUT'
            output_param = 'OUTPUT'
        else:
            input_param = 'input'
            output_param = 'output'
        if output:
            base_cmd += '--{output_param} {output} '.format(output_param=output_param, output=output)
        if input_bam:
            base_cmd += '--{input_param} {input_bam} '.format(input_param=input_param, input_bam=input_bam)

        if spark_core > 1 or 'Spark' in run_cls:
            base_cmd += '--spark-master local[{spark_core}] ' \
                        '--conf spark.driver.host=localhost ' \
                        '--conf spark.network.timeout=800 ' \
                        '--conf spark.local.dir={tmp_dir} '
        if reference == 'default':
            reference = self.dataset.reference_genome
        if reference:
            base_cmd += ' --reference ' + reference
        if ext:
            base_cmd += ' ' + ext
        return base_cmd.format(
            tmp_dir=self.gatk_run_dir,
            spark_core=spark_core,
            memory=memory,
            run_cls=run_cls,
            output=output
        )

    def gatk_cmd(self, run_cls, output, memory=2, input=None, reference='default', spark_core=1, ext=None):
        return self._gatk_cmd(run_cls, output=output, memory=memory, input=input, reference=reference,
                              spark_core=spark_core, ext=ext)

    def gatk_picard_cmd(self, run_cls, output, memory=2, input=None, reference='default', ext=None):
        return self._gatk_cmd(run_cls, output=output, memory=memory, input=input, reference=reference,
                              ext=ext, picard_tool=True)

    @property
    def basename(self):
        return os.path.join(self.gatk_run_dir, self.dataset.user_sample_id)

    @property
    def sorted_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')

    @property
    def recal_bam(self):
        return self.basename + '_recal.bam'

    @property
    def output_grp(self):
        return self.basename + '.grp'

    @property
    def output_intervals(self):
        return self.basename + '.intervals'

    @property
    def sample_gvcf(self):
        return self.basename + '.g.vcf.gz'

    @property
    def genotyped_vcf(self):
        return self.basename + '.vcf.gz'

    @property
    def raw_snp_vcf(self):
        return self.basename + '_raw_snp.vcf'

    @property
    def filter_snp_vcf(self):
        return self.basename + '_filter_snp.vcf'

    @property
    def dbsnp(self):
        dbsnp = cfg.query('genomes', self.dataset.genome_version, 'dbsnp')
        if not dbsnp and self.dataset.pipeline.name == 'variant_calling':
            raise AnalysisDriverError('Could not find dbsnp file for %s' % self.dataset.name)
        return dbsnp

    @property
    def known_indels(self):
        return cfg.query('genomes', self.dataset.genome_version, 'known_indels')


class PostAlignmentScatter(GATK4Stage):

    @property
    def split_file_dir(self):
        return os.path.join(self.job_dir, 'post_alignment_split')

    def split_genome_files(self):
        """
        This function create the bed file representing all the chunks of genomes required.
        It also create a dictionary where keys are the first chunk of each bed file and values are the file path.
        """
        os.makedirs(self.split_file_dir, exist_ok=True)
        chunk_to_file = {}
        for chunks in self.split_genome_in_chunks():
            chunk_file = os.path.join(self.split_file_dir, self.dataset.name + '_region_%s-%s-%s.bed' % chunks[0])
            if not os.path.exists(chunk_file):
                with open(chunk_file, 'w') as open_file:
                    for chunk in chunks:
                        open_file.write('\t'.join([str(e) for e in chunk]) + '\n')
            chunk_to_file[chunks[0]] = chunk_file
        return chunk_to_file

    def split_genome_in_chunks(self):
        """
        Split a genome in chunks of a specific size or at the end of chromosome and the aggregate the small chunks to
        roughly match the other chunk size.
        """
        fai_file = self.dataset.reference_genome + '.fai'
        chunk_size = 20000000
        with open(fai_file) as open_file:
            chunks = []
            for line in open_file:
                sp_line = line.strip().split()
                chromosome = sp_line[0]
                length = int(sp_line[1])
                last = 0
                for chunki in range(length // chunk_size + min(1, length % chunk_size)):
                    new = last + chunk_size
                    chunks.append((chromosome, last, min(new, length)))
                    last = new
        # merge small consecutive chunks together
        final_chunks = []
        current_chunks = []
        final_chunks.append(current_chunks)
        current_chunks_size = 0
        for chunk in chunks:
            if current_chunks_size + (chunk[2] - chunk[1]) > chunk_size:
                current_chunks = []
                current_chunks_size = 0
                final_chunks.append(current_chunks)
            current_chunks.append(chunk)
            current_chunks_size += chunk[2] - chunk[1]
        return final_chunks

    def split_genome_chromosomes(self, with_unmapped=False):
        """
        Split the genome per chromosomes aggregating smaller chromosomes to similar length as the longuest chromsome
        Code inspired from GATK best practice workflow:
        https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/blob/190945e358a6ee7a8c65bacd7b28c66527383376/PairedEndSingleSampleWf.wdl#L969

        :return: list of list of chromosome names
        """
        fai_file = self.dataset.reference_genome + '.fai'
        with open(fai_file) as open_file:
            sequence_tuple_list = []
            for line in open_file:
                sp_line = line.strip().split()
                sequence_tuple_list.append((sp_line[0], int(sp_line[1])))
            longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
        chunks = []
        current_chunks = []
        chunks.append(current_chunks)
        temp_size = 0
        for sequence_tuple in sequence_tuple_list:
            if temp_size + sequence_tuple[1] <= longest_sequence:
                temp_size += sequence_tuple[1]
                current_chunks.append(sequence_tuple[0])
            else:
                current_chunks = []
                chunks.append(current_chunks)
                current_chunks.append(sequence_tuple[0])
                temp_size = sequence_tuple[1]
        # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
        if with_unmapped:
            chunks.append(['unmapped'])
        return chunks

    def split_base_recal_grp(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_base_recal_grp_%s.grp' % chunk)

    def split_recal_bam(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_recal_%s.bam' % chunk)

    def gvcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_haplotype_caller_%s-%s-%s.g.vcf.gz' % chunk)

    def vcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_genotype_gvcf_%s-%s-%s.vcf.gz' % chunk)


class ScatterBaseRecalibrator(PostAlignmentScatter):

    def base_recalibrator_cmd(self, chrom_names):
        return self.gatk_cmd(
            'BaseRecalibrator', self.split_base_recal_grp(chrom_names[0]), input=self.sorted_bam,
            memory=6, ext='--known-sites  ' + self.dbsnp + ' --intervals ' + ' --intervals '.join(chrom_names)
        )

    def _run(self):
        return executor.execute(
            *[self.base_recalibrator_cmd(chrom_names)
              for chrom_names in self.split_genome_chromosomes()],
            job_name='gatk_base_recal',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherBQSRReport(PostAlignmentScatter):
    def _run(self):
        bqsr_reports_list = os.path.join(self.split_file_dir, self.dataset.name + 'bqsr_reports.list')
        with open(bqsr_reports_list, 'w') as open_file:
            for chrom_names in self.split_genome_chromosomes():
                open_file.write(self.split_base_recal_grp(chrom_names[0]) + '\n')

        gather_bqsr_status = executor.execute(
            self.gatk_cmd('GatherBQSRReports', self.output_grp, input=bqsr_reports_list, memory=6, reference=None),
            job_name='gather_bqsr',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return gather_bqsr_status


class ScatterApplyBQSR(PostAlignmentScatter):

    def apply_bqsr_cmd(self, chrom_names):
        return self.gatk_cmd(
            'ApplyBQSR', self.split_recal_bam(chrom_names[0]), input=self.sorted_bam,
            memory=6, ext='--bqsr-recal-file ' + self.output_grp + ' --jdk-deflater --jdk-inflater '
                          '--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 '
                          '--static-quantized-quals 40' + ' --intervals ' + ' --intervals '.join(chrom_names)
        )

    def _run(self):
        return executor.execute(
            *[self.apply_bqsr_cmd(chrom_names)
              for chrom_names in self.split_genome_chromosomes(with_unmapped=True)],
            job_name='apply_bqsr',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherRecalBam(PostAlignmentScatter):

    def _run(self):
        bam_file_list = os.path.join(self.split_file_dir, self.dataset.name + 'recal_bam.list')
        with open(bam_file_list, 'w') as open_file:
            for chrom_names in self.split_genome_chromosomes(with_unmapped=True):
                open_file.write(self.split_recal_bam(chrom_names[0]) + '\n')

        gather_bam_status = executor.execute(
            picard_command('GatherBamFiles', input_file=bam_file_list, output_file=self.recal_bam,
                           tmp_dir=self.gatk_run_dir, memory=6, assume_sorted=False,
                           picard_params={'CREATE_INDEX': 'true'}),
            job_name='gather_recal_bam',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return gather_bam_status


class SplitHaplotypeCaller(PostAlignmentScatter):

    def haplotype_caller_cmd(self, chunk, region_file):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.gvcf_per_chunk(chunk),
            input=self.recal_bam,
            memory=6,
            spark_core=1,
            ext=' --sample-ploidy 2 --emit-ref-confidence GVCF --intervals ' + region_file
        )
        for annot in ('BaseQualityRankSumTest', 'ClippingRankSumTest', 'Coverage', 'DepthPerAlleleBySample',
                      'DepthPerSampleHC', 'FisherStrand', 'MappingQuality', 'MappingQualityRankSumTest',
                      'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest', 'RMSMappingQuality'):
            haplotype_cmd += ' --annotation ' + annot
        if self.dbsnp:
            haplotype_cmd += ' --dbsnp ' + self.dbsnp
        if self.dataset.library_type == 'pcr-free':
            haplotype_cmd += '--pcr-indel-model NONE '

        return haplotype_cmd

    def _run(self):
        haplotype_status = executor.execute(
            *[self.haplotype_caller_cmd(chunk, region_file)
              for chunk, region_file in self.split_genome_files().items()],
            job_name='split_hap_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return haplotype_status


class GatherGVCF(PostAlignmentScatter):
    def _run(self):
        gvcf_list = os.path.join(self.split_file_dir, self.dataset.name + '_g.vcf.list')
        with open(gvcf_list, 'w') as open_file:
            for chunks in self.split_genome_in_chunks():
                open_file.write(self.gvcf_per_chunk(chunks[0]) + '\n')

        concat_vcf_status = executor.execute(
            self.gatk_picard_cmd('GatherVcfs', self.sample_gvcf, input=gvcf_list, memory=6, reference=None),
            job_name='gather_hap_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        if concat_vcf_status == 0:
            concat_vcf_status = tabix_vcf(self.exec_dir, self.genotyped_vcf)

        return concat_vcf_status


class SplitGenotypeGVCFs(PostAlignmentScatter):

    def genotypegvcf_cmd(self, chunk, region_file):
        genotypegvcf_cmd = self.gatk_cmd(
            'GenotypeGVCFs',
            self.vcf_per_chunk(chunk),
            memory=6,
            spark_core=1,
            ext='--variant ' + self.gvcf_per_chunk(chunk) + ' --sample-ploidy 2 --intervals ' + region_file
        )
        for annot in ('BaseQualityRankSumTest', 'ClippingRankSumTest', 'Coverage', 'DepthPerAlleleBySample',
                      'DepthPerSampleHC', 'FisherStrand', 'MappingQuality', 'MappingQualityRankSumTest',
                      'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest', 'RMSMappingQuality'):
            genotypegvcf_cmd += ' --annotation ' + annot

        return genotypegvcf_cmd

    def _run(self):
        return executor.execute(
            *[self.genotypegvcf_cmd(chunk, region_file)
              for chunk, region_file in self.split_genome_files().items()],
            job_name='split_genotype_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()


class GatherVCF(PostAlignmentScatter):

    def _run(self):
        vcf_list = os.path.join(self.split_file_dir, self.dataset.name + '.vcf.list')
        with open(vcf_list, 'w') as open_file:
            for chunks in self.split_genome_in_chunks():
                open_file.write(self.vcf_per_chunk(chunks[0]) + '\n')

        concat_vcf_status = executor.execute(
            self.gatk_picard_cmd('GatherVcfs', self.genotyped_vcf, input=vcf_list, memory=6, reference=None),
            job_name='gather_geno_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        if concat_vcf_status == 0:
            concat_vcf_status = tabix_vcf(self.exec_dir, self.genotyped_vcf)

        return concat_vcf_status


class GATK4VariantCall(GATK4Stage):
    @property
    def vqsr_datasets(self):
        vqsr_datasets = cfg.query('genomes', self.dataset.genome_version, 'vqsr')
        if not vqsr_datasets:
            raise AnalysisDriverError(
                'Could not find VQSR training and evaluation sets file for genome %s' % self.dataset.genome_version)
        return vqsr_datasets

    @property
    def vqsr_snps_output_recall(self):
        return self.basename + '_vqsr_snps_output.recall'

    @property
    def vqsr_snps_tranches(self):
        return self.basename + '_vqsr_snps.tranches'

    @property
    def vqsr_snps_R_script(self):
        return self.basename + '_vqsr_snps.R'

    @property
    def vqsr_filtered_snps_vcf(self):
        return self.basename + '_vqsr_snps.vcf.gz'

    @property
    def vqsr_indels_output_recall(self):
        return self.basename + '_vqsr_indels_output.recall'

    @property
    def vqsr_indel_tranches(self):
        return self.basename + '_vqsr_indels_tranches'

    @property
    def vqsr_indel_R_script(self):
        return self.basename + '_vqsr_indels.R'

    @property
    def vqsr_filtered_indel_vcf(self):
        return self.basename + '_vqsr_indels.vcf.gz'

    @property
    def vqsr_filtered_vcf(self):
        return self.basename + '_vqsr.vcf.gz'


class VQSRFiltrationSNPs(GATK4VariantCall):
    def _run(self):
        vqsr_datasets = self.vqsr_datasets
        cmd = self.gatk_cmd(
            'VariantRecalibrator',
            self.vqsr_snps_output_recall,
            memory=16,
            ext='-V {input_vcf} '
                '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} '
                '--resource:omni,known=false,training=true,truth=false,prior=12.0 {omni} '
                '--resource:1000G,known=false,training=true,truth=false,prior=10.0 {thousand_genomes} '
                '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} '
                '-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '
                '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP '
                '--tranches-file {ouput_tranches} --rscript-file {output_R_script}'.format(
                    input_vcf=self.genotyped_vcf,
                    hapmap=vqsr_datasets.get('hapmap'),
                    omni=vqsr_datasets.get('omni'),
                    thousand_genomes=vqsr_datasets.get('thousand_genomes'),
                    dbsnp=vqsr_datasets.get('dbsnp'),
                    ouput_tranches=self.vqsr_snps_tranches,
                    output_R_script=self.vqsr_snps_R_script
            )
        )
        return executor.execute(
            cmd,
            job_name='vqsr_snp',
            working_dir=self.exec_dir,
            cpu=1,
            mem=16
        ).join()


class VQSRFiltrationIndels(GATK4VariantCall):
    def _run(self):
        vqsr_datasets = self.vqsr_datasets
        cmd = self.gatk_cmd(
            'VariantRecalibrator',
            self.vqsr_indels_output_recall,
            memory=16,
            ext='-V {input_vcf} '
                '--resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} '
                '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} '
                '-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 '
                '-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode INDEL '
                '--tranches-file {ouput_tranches} --rscript-file {output_R_script}'.format(
                    input_vcf=self.genotyped_vcf,
                    mills=vqsr_datasets.get('mills'),
                    dbsnp=vqsr_datasets.get('dbsnp'),
                    ouput_tranches=self.vqsr_indel_tranches,
                    output_R_script=self.vqsr_indel_R_script,
                )
        )
        return executor.execute(
            cmd,
            job_name='vqsr_indel',
            working_dir=self.exec_dir,
            cpu=1,
            mem=16
        ).join()


class ApplyVQSRSNPs(GATK4VariantCall):
    def _run(self):
        cmd = self.gatk_cmd(
            'ApplyVQSR',
            self.vqsr_filtered_snps_vcf,
            memory=16,
            ext='-V {input_vcf} -mode SNP --tranches-file {ouput_tranches} --truth-sensitivity-filter-level 99.0 '
                '--recal-file {recal_file}'.format(
                input_vcf=self.genotyped_vcf,
                ouput_tranches=self.vqsr_snps_tranches,
                recal_file=self.vqsr_snps_output_recall
            )
        )
        return executor.execute(
            cmd,
            job_name='apply_vqsr_snps',
            working_dir=self.exec_dir,
            cpu=1,
            mem=16
        ).join()


class ApplyVQSRIndels(GATK4VariantCall):
    def _run(self):
        cmd = self.gatk_cmd(
            'ApplyVQSR',
            self.vqsr_filtered_indels_vcf,
            memory=8,
            ext='-V {input_vcf} -mode INDELS --tranches-file {ouput_tranches} --truth-sensitivity-filter-level 99.0 '
                '--recal-file {recal_file}'.format(
                input_vcf=self.genotyped_vcf,
                ouput_tranches=self.vqsr_indels_tranches,
                recal_file=self.vqsr_indels_output_recall
            )
        )
        return executor.execute(
            cmd,
            job_name='apply_vqsr_indels',
            working_dir=self.exec_dir,
            cpu=1,
            mem=16
        ).join()


class MergeVariants(GATK4VariantCall):

    def _run(self):
        vcf_list = os.path.join(self.split_file_dir, self.dataset.name + 'vqsr_filtered.vcf.list')
        with open(vcf_list, 'w') as open_file:
            open_file.write(self.vqsr_filtered_snps_vcf + '\n')
            open_file.write(self.vqsr_filtered_indels_vcf + '\n')
        cmd = self.gatk_picard_cmd('MergeVcfs', self.vqsr_filtered_vcf, input=vcf_list)
        return executor.execute(
            cmd,
            job_name='merge_vqsr',
            working_dir=self.exec_dir,
            cpu=1,
            mem=8
        ).join()


class VariantAnnotation(GATK4VariantCall):
    """Annotate a vcf file using snpEff."""

    input_vcf = Parameter(default=None)

    @property
    def snpeff_basename(self):
        base, ext = os.path.splitext(self.input_vcf)
        if ext == '.gz':
            base, ext = os.path.splitext(base)
        return base

    @property
    def snps_effects_csv(self):
        return self.snpeff_basename + '_snpseff.csv'

    @property
    def snps_effects_html(self):
        return self.snpeff_basename + '_snpseff.html'

    @property
    def output_vcf(self):
        return self.snpeff_basename + '_effects' + '.vcf.gz'

    def _run(self):
        cmd = ('{snpEff} -Xmx{memory}g -Djava.io.tmpdir={tmp_dir} eff '
               '-dataDir {snpEff_ressource} -hgvs -noLog -i vcf -o vcf '
               '-csvStats {effects_csv} -s {effects_html} {database_name} {input_vcf)'
               ' | {bgzip} --threads 16 -c > {output_vcf}').format(
            snpEff=toolset['snpEff'],
            bgzip=toolset['bgzip'],
            memory=20,
            tmp_dir=self.job_dir,
            snpEff_ressource=cfg.query('snpEff','ressource'),
            effects_csv=self.snps_effects_csv,
            effects_html=self.snps_effects_html,
            database_name=cfg.query('genomes', self.dataset.genome_version, 'snpEff'),
            input_vcf=self.input_vcf,
            output_vcf=self.output_vcf
        )
        return executor.execute(
            cmd,
            job_name='snpeff',
            working_dir=self.exec_dir,
            mem=20
        ).join()


class SelectVariants(GATK4Stage):
    def _run(self):
        select_var_command = self.gatk_cmd('SelectVariants', self.raw_snp_vcf, nct=1, nt=16)
        select_var_command += ' -V ' + self.genotyped_vcf + '.gz'
        select_var_command += ' -selectType SNP '
        select_variants_status = executor.exec_dir(
            select_var_command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()
        return select_variants_status


class SNPsFiltration(GATK4Stage):
    def _run(self):
        filter_array = [
            'QD < 2.0',
            'FS > 60.0',
            'MQ < 40.0',
            'MQRankSum < -12.5',
            'ReadPosRankSum < -8.0',
            'SOR > 3.0'
        ]
        filters = "'" + ' || '.join(filter_array) + "'"
        var_filter_command = self.gatk_cmd('VariantFiltration', self.filter_snp_vcf, nct=1, nt=1)
        var_filter_command += " -V " + self.raw_snp_vcf + '.gz'
        var_filter_command += " --filterExpression " + filters
        var_filter_command += " --filterName 'SNP_HARD_FILTER'"
        variant_filter_status = executor.exec_dir(
            var_filter_command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()
        return variant_filter_status


class IndelsFiltration(GATK4Stage):
    def _run(self):
        filter_array = [
            'QD < 2.0',
            'FS > 200.0',
            'ReadPosRankSum < -20.0',
            'SOR > 10.0'
        ]
        filters = "'" + ' || '.join(filter_array) + "'"
        var_filter_command = self.gatk_cmd('VariantFiltration', self.filter_snp_vcf, nct=1, nt=1)
        var_filter_command += " -V " + self.raw_snp_vcf + '.gz'
        var_filter_command += " --filterExpression " + filters
        var_filter_command += " --filterName 'INDEL_HARD_FILTER'"
        variant_filter_status = executor.exec_dir(
            var_filter_command,
            job_name='var_filtration',
            working_dir=self.gatk_run_dir,
            mem=16
        ).join()
        return variant_filter_status



# class BamFileGenerationPipeline:
#
#     toolset_type = 'gatk4_sample_processing'
#     name = 'variant_calling_gatk4'
#
#     def __init__(self, dataset):
#         self.dataset =dataset
#
#     def stage(self, klass, **params):
#         return klass(dataset=self.dataset, **params)
#
#     def pipeline(self):
#         fastq_index = self.stage(FastqIndex)
#         split_bwa = self.stage(SplitBWA, previous_stages=[fastq_index])
#         merge_bam = self.stage(MergeBam, previous_stages=[split_bwa])
#         return merge_bam
#
#
# class GATK4VariantCall(BamFileGenerationPipeline):
#
#     def pipeline(self):
#         merge_bam = super().pipeline()
#         contam = self.stage(qc.FastqScreen, previous_stages=[merge_bam], fq_pattern=bwa.fq_pattern)
#         blast = self.stage(qc.Blast, previous_stages=[bwa], fastq_file=bwa.fq_pattern.replace('?', '1'))
#         samtools_stat = self.stage(SamtoolsStats, previous_stages=[fastqc, bwa, contam, blast])
#         samtools_depth = self.stage(qc.SamtoolsDepth, bam_file=bwa.exp_bam_path, previous_stages=[bwa])
#         base_recal = self.stage(BaseRecal, previous_stages=[merge_bam])
#         apply_bqsr = self.stage(ApplyBQSR, previous_stages=[base_recal])
#         return apply_bqsr


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    # merge fastq to do contamination check
    merge_fastqs = stage(MergeFastqs)
    contam = stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
    blast = stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))
    geno_val = stage(qc.GenotypeValidation, fq_pattern=merge_fastqs.fq_pattern, previous_stages=[merge_fastqs])

    # create fastq index then align and recalibrate via scatter-gather strategy
    fastq_index = stage(FastqIndex)
    split_bwa = stage(SplitBWA, previous_stages=[fastq_index])
    merge_bam_dup = stage(MergeBamAndDup, previous_stages=[split_bwa])
    base_recal = stage(ScatterBaseRecalibrator, previous_stages=[merge_bam_dup])
    gather_bqsr = stage(GatherBQSRReport, previous_stages=[base_recal])
    apply_bqsr = stage(ScatterApplyBQSR, previous_stages=[gather_bqsr])
    merge_bam = stage(GatherRecalBam, previous_stages=[apply_bqsr])

    # bam file QC
    verify_bam_id = stage(qc.VerifyBamID, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])
    samtools_stat = stage(SamtoolsStats, previous_stages=[merge_bam])
    samtools_depth = stage(qc.SamtoolsDepth, bam_file=merge_bam.recal_bam, previous_stages=[merge_bam])

    # variants call via scatter-gather strategy
    haplotype_caller = stage(SplitHaplotypeCaller, previous_stages=[merge_bam])
    gather_gcvf = stage(GatherGVCF, previous_stages=[haplotype_caller])
    genotype_gcvf = stage(SplitGenotypeGVCFs, previous_stages=[haplotype_caller])
    gather_vcf = stage(GatherVCF, previous_stages=[genotype_gcvf])

    # variant filtering
    filter_snps = stage(VQSRFiltrationSNPs, previous_stages=[gather_vcf])
    apply_vqsr_snps = stage(ApplyVQSRSNPs, previous_stages=[filter_snps])
    filter_indels = stage(VQSRFiltrationIndels, previous_stages=[gather_vcf])
    apply_vqsr_indels = stage(ApplyVQSRIndels, previous_stages=[filter_indels])

    # final variant merge
    merge_vqsr = stage(MergeVariants, previous_stages=[apply_vqsr_snps, apply_vqsr_indels])

    # variant file QC
    gender_val = stage(qc.GenderValidation, vcf_file=merge_vqsr.vqsr_filtered_vcf, previous_stages=[merge_vqsr])
    vcfstats = stage(qc.VCFStats, vcf_file=merge_vqsr.vqsr_filtered_vcf, previous_stages=[merge_vqsr])

    final_stages = [contam, blast, geno_val, gender_val, vcfstats, verify_bam_id, samtools_depth, samtools_stat,
                    gather_gcvf]

    # output = stage(common.SampleDataOutput, previous_stages=post_alignment_qc, output_fileset='bcbio')
    # cleanup = stage(common.Cleanup, previous_stages=[output])
    # review = stage(common.SampleReview, previous_stages=[cleanup])

    return final_stages

