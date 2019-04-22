import os
import shutil
from os.path import join

from egcg_core import executor, util
from egcg_core.config import cfg
from egcg_core.constants import ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_ELEMENT_ID, ELEMENT_PROJECT_ID, \
    ELEMENT_RUN_NAME, ELEMENT_LANE

from analysis_driver import segmentation
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver import quality_control as qc
from analysis_driver.pipelines import common
from analysis_driver.pipelines.common import tabix_vcf, MergeFastqs, SamtoolsStats
from analysis_driver.segmentation import Parameter, ListParameter
from analysis_driver.tool_versioning import toolset

toolset_type = 'gatk4_sample_processing'
name = 'qc_gatk4'


class GATK4FilePath(segmentation.Stage):

    @property
    def gatk_run_dir(self):
        d = os.path.join(self.job_dir, 'gatk4')
        os.makedirs(d, exist_ok=True)
        return d

    @property
    def gatk4_basename(self):
        return os.path.join(self.gatk_run_dir, self.dataset.user_sample_id)

    @property
    def sorted_bam(self):
        return os.path.join(self.job_dir, self.dataset.name + '.bam')

    @property
    def recal_bam(self):
        return self.gatk4_basename + '_recal.bam'

    @property
    def output_grp(self):
        return self.gatk4_basename + '.grp'

    @property
    def output_intervals(self):
        return self.gatk4_basename + '.intervals'

    @property
    def sample_gvcf(self):
        return self.gatk4_basename + '.g.vcf.gz'


    @property
    def raw_snps_vcf(self):
        return self.gatk4_basename + '_raw_snp.vcf'

    @property
    def raw_indels_vcf(self):
        return self.gatk4_basename + '_raw_indel.vcf'

    @property
    def hard_filtered_snps_vcf(self):
        return self.gatk4_basename + '_hard_snps.vcf.gz'

    @property
    def hard_filtered_indels_vcf(self):
        return self.gatk4_basename + '_hard_filter_indels.vcf.gz'

    @property
    def hard_filtered_vcf(self):
        return self.gatk4_basename + '_hard_filter.vcf.gz'

    @property
    def genotyped_vcf(self):
        return self.gatk4_basename + '.vcf.gz'

    @property
    def dbsnp(self):
        return cfg.query('genomes', self.dataset.genome_version, 'dbsnp')

    @property
    def snps_effects_csv(self):
        return self.gatk4_basename + '_snpseff.csv'

    @property
    def snps_effects_html(self):
        return self.gatk4_basename + '_snpseff.html'

    @property
    def snps_effects_output_vcf(self):
        return self.gatk4_basename + '_snpseff' + '.vcf.gz'

    @property
    def vqsr_datasets(self):
        vqsr_datasets = cfg.query('genomes', self.dataset.genome_version, 'vqsr')
        if not vqsr_datasets:
            raise AnalysisDriverError(
                'Could not find VQSR training and evaluation sets file for genome %s' % self.dataset.genome_version)
        return vqsr_datasets

    @property
    def vqsr_snps_output_recall(self):
        return self.gatk4_basename + '_vqsr_snps_recall.vcf.gz'

    @property
    def vqsr_snps_tranches(self):
        return self.gatk4_basename + '_vqsr_snps.tranches'

    @property
    def vqsr_snps_R_script(self):
        return self.gatk4_basename + '_vqsr_snps.R'

    @property
    def vqsr_filtered_snps_vcf(self):
        return self.gatk4_basename + '_vqsr_snps.vcf.gz'

    @property
    def vqsr_indels_output_recall(self):
        return self.gatk4_basename + '_vqsr_indels_recall.vcf.gz'

    @property
    def vqsr_indels_tranches(self):
        return self.gatk4_basename + '_vqsr_indels_tranches'

    @property
    def vqsr_indels_R_script(self):
        return self.gatk4_basename + '_vqsr_indels.R'

    @property
    def vqsr_filtered_indels_vcf(self):
        return self.gatk4_basename + '_vqsr_indels.vcf.gz'

    @property
    def vqsr_filtered_vcf(self):
        return self.gatk4_basename + '_vqsr.vcf.gz'


class SplitFastqStage(GATK4FilePath):
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
        cmd_template = 'gunzip -c {ipath} | {pbgzip} -n 16  -c /dev/stdin > {opath}'
        for run_element in self.dataset.run_elements:
            if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                for ifastqs, ofastqs in zip(
                        self._find_fastqs_for_run_element(run_element),
                        self._indexed_fastq_for_run_element(run_element)
                ):
                    commands.append(cmd_template.format(ipath=ifastqs, pbgzip=toolset['pbgzip'], opath=ofastqs))
        exit_status = executor.execute(
            *commands,
            job_name='compress_fastq',
            working_dir=self.exec_dir,
            mem=8,
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
                working_dir=self.exec_dir
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
            working_dir=self.exec_dir,
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

    def merge_command(self):
        cat_cmd = '{bamcat} level=0 tmpfile={cat_tmp} `cat {bam_list_file}`'.format(
            bamcat=toolset['bamcat'],
            cat_tmp=os.path.join(self.create_tmp_dir(), self.dataset.name),
            bam_list_file=self.all_chunk_bam_list_file()
        )
        bamsormadup_cmd = '{bamsormadup} threads=16 tmpfile={dup_tmp} indexfilename={merged_bam}.bai > {merged_bam}'.format(
            bamsormadup=toolset['biobambam_sortmapdup'],
            dup_tmp=os.path.join(self.create_tmp_dir(), self.dataset.name),
            merged_bam=self.sorted_bam
        )
        return 'set -o pipefail; ' + ' | '.join([cat_cmd, bamsormadup_cmd])

    def _run(self):
        exit_status = executor.execute(
            self.merge_command(),
            job_name='merge_dup_bam',
            working_dir=self.exec_dir,
            cpus=6,
            mem=36
        ).join()
        return exit_status


class GATK4Stage(GATK4FilePath):

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
        if input:
            base_cmd += '--{input_param} {input} '.format(input_param=input_param, input=input)

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
            tmp_dir=self.create_tmp_dir(),
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


class PostAlignmentScatter(GATK4Stage):

    @property
    def split_file_dir(self):
        d = os.path.join(self.job_dir, 'post_alignment_split')
        os.makedirs(d, exist_ok=True)
        return d

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
        max_nb_contig_per_chunk = 1000
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
        # merge small consecutive chunks together with an upper limit set by max_nb_contig_per_chunk
        final_chunks = []
        current_chunks = []
        final_chunks.append(current_chunks)
        current_chunks_size = 0
        for chunk in chunks:
            if current_chunks_size + (chunk[2] - chunk[1]) > chunk_size or len(current_chunks) > max_nb_contig_per_chunk:
                current_chunks = []
                current_chunks_size = 0
                final_chunks.append(current_chunks)
            current_chunks.append(chunk)
            current_chunks_size += chunk[2] - chunk[1]
        return final_chunks

    def vcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_haplotype_caller_%s-%s-%s.vcf.gz' % chunk)


class SplitHaplotypeCaller(PostAlignmentScatter):

    bam_file = Parameter()

    def haplotype_caller_cmd(self, chunk, region_file):
        haplotype_cmd = self.gatk_cmd(
            'HaplotypeCaller',
            self.vcf_per_chunk(chunk),
            input=self.bam_file,
            memory=6,
            spark_core=1,
            ext=' --sample-ploidy 2 --intervals ' + region_file
        )
        for annot in ('BaseQualityRankSumTest', 'ClippingRankSumTest', 'Coverage', 'DepthPerAlleleBySample',
                      'DepthPerSampleHC', 'FisherStrand', 'MappingQuality', 'MappingQualityRankSumTest',
                      'MappingQualityZero', 'QualByDepth', 'ReadPosRankSumTest', 'RMSMappingQuality'):
            haplotype_cmd += ' --annotation ' + annot
        if self.dbsnp:
            haplotype_cmd += ' --dbsnp ' + self.dbsnp
        if self.dataset.library_preparation == 'pcr-free':
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


class GatherVCF(PostAlignmentScatter):

    def _run(self):
        vcf_list = os.path.join(self.split_file_dir, self.dataset.name + '_vcf.list')
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


class SelectSNPs(GATK4Stage):

    input_vcf = Parameter()

    def _run(self):
        select_var_command = self.gatk_cmd('SelectVariants', self.raw_snps_vcf)
        select_var_command += ' -V ' + self.input_vcf
        select_var_command += ' -select-type SNP '
        select_variants_status = executor.execute(
            select_var_command,
            job_name='snp_select',
            working_dir=self.exec_dir,
            mem=16
        ).join()
        return select_variants_status


class SelectIndels(GATK4Stage):

    input_vcf = Parameter()

    def _run(self):
        select_var_command = self.gatk_cmd('SelectVariants', self.raw_indels_vcf)
        select_var_command += ' -V ' + self.input_vcf
        select_var_command += ' -select-type INDEL '
        select_variants_status = executor.execute(
            select_var_command,
            job_name='indel_select',
            working_dir=self.exec_dir,
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
        var_filter_command = self.gatk_cmd('VariantFiltration', self.hard_filtered_snps_vcf)
        var_filter_command += " -V " + self.raw_snps_vcf
        var_filter_command += " --filter-expression " + filters
        var_filter_command += " --filter-name 'SNP_HARD_FILTER'"
        variant_filter_status = executor.execute(
            var_filter_command,
            job_name='snps_filtration',
            working_dir=self.exec_dir,
            mem=8
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
        var_filter_command = self.gatk_cmd('VariantFiltration', self.hard_filtered_indels_vcf)
        var_filter_command += " -V " + self.raw_indels_vcf
        var_filter_command += " --filter-expression " + filters
        var_filter_command += " --filter-name 'INDEL_HARD_FILTER'"
        variant_filter_status = executor.execute(
            var_filter_command,
            job_name='indel_filtration',
            working_dir=self.exec_dir,
            mem=16
        ).join()
        return variant_filter_status


class MergeVariants(GATK4Stage):
    vcf_files = ListParameter()
    output_vcf_file = Parameter()

    def _run(self):
        # get the name of the list from the name of the first vcf to avoid overwriting another file
        base_name = os.path.splitext(os.path.basename(self.vcf_files[0]))[0]
        vcf_list = os.path.join(self.job_dir, base_name + '.list')
        with open(vcf_list, 'w') as open_file:
            for f in self.vcf_files:
                open_file.write(f + '\n')
        cmd = self.gatk_picard_cmd('MergeVcfs', self.output_vcf_file, input=vcf_list, reference=None)
        return executor.execute(
            cmd,
            job_name='merge_vqsr',
            working_dir=self.exec_dir,
            cpus=1,
            mem=8
        ).join()

def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    # merge fastq to do contamination check
    merge_fastqs = stage(MergeFastqs)
    contam = stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
    blast = stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))

    # create fastq index then align via scatter-gather strategy
    fastq_index = stage(FastqIndex)
    split_bwa = stage(SplitBWA, previous_stages=[fastq_index])
    merge_bam_dup = stage(MergeBamAndDup, previous_stages=[split_bwa])

    # bam file QC
    samtools_stat = stage(SamtoolsStats, bam_file=merge_bam_dup.sorted_bam, previous_stages=[merge_bam_dup])
    samtools_depth = stage(qc.SamtoolsDepth, bam_file=merge_bam_dup.sorted_bam, previous_stages=[merge_bam_dup])

    # variants call via scatter-gather strategy
    haplotype_caller = stage(SplitHaplotypeCaller, bam_file=merge_bam_dup.sorted_bam, previous_stages=[merge_bam_dup])
    gather_vcf = stage(GatherVCF, previous_stages=[haplotype_caller])

    # variant filtering
    select_snps = stage(SelectSNPs, input_vcf=gather_vcf.genotyped_vcf, previous_stages=[gather_vcf])
    select_indels = stage(SelectIndels, input_vcf=gather_vcf.genotyped_vcf, previous_stages=[gather_vcf])
    filter_snps = stage(SNPsFiltration, previous_stages=[select_snps])
    filter_indels = stage(IndelsFiltration, previous_stages=[select_indels])

    # final variant merge
    merge_vqsr = stage(MergeVariants,
                       vcf_files=[filter_snps.hard_filtered_snps_vcf, filter_indels.hard_filtered_indels_vcf],
                       output_vcf_file=filter_indels.hard_filtered_vcf,
                       previous_stages=[filter_snps, filter_indels])

    # variant file QC
    vcfstats = stage(qc.VCFStats, vcf_file=merge_vqsr.hard_filtered_vcf, previous_stages=[merge_vqsr])

    final_stages = [contam, blast, vcfstats, samtools_depth, samtools_stat]

    output = stage(common.SampleDataOutput, previous_stages=final_stages, output_fileset='gatk4_qc')
    review = stage(common.SampleReview, previous_stages=[output])
    # cleanup = stage(common.Cleanup, previous_stages=[output])

    return [review]

