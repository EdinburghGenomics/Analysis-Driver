import os
from egcg_core import executor, util
from egcg_core.config import cfg
from egcg_core.constants import ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_ELEMENT_ID, ELEMENT_PROJECT_ID, \
    ELEMENT_RUN_NAME, ELEMENT_LANE

from analysis_driver import quality_control as qc
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.pipelines import Pipeline, common
from analysis_driver.segmentation import Parameter, ListParameter, Stage
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.helper_functions import split_in_chunks


class ChunkHandler:
    def __init__(self, dataset):
        self.dataset = dataset
        self.genome_chunks = self.split_genome_in_chunks()

    def split_genome_in_chunks(self):
        """
        Split a genome in chunks of a specific size or at the end of chromosome and then aggregate the small chunks to
        roughly match the other chunk size.
        :return: list of list
        """
        fai_file = self.dataset.reference_genome + '.fai'
        chunk_size = 20000000
        max_nb_contig_per_chunk = 1000
        with open(fai_file) as open_file:
            chunks = []
            for line in open_file:
                sp_line = line.strip().split()
                length = int(sp_line[1])
                # The indices in bed files are 0 based and end excluded
                # https://genome.ucsc.edu/FAQ/FAQformat.html#format1
                chunks.extend([
                    (sp_line[0], start, end)
                    for start, end in split_in_chunks(length, chunk_size, zero_based=True, end_inclusive=False)
                ])
        # merge small consecutive chunks together with an upper limit set by max_nb_contig_per_chunk
        final_chunks = []
        current_chunks = []
        final_chunks.append(current_chunks)
        current_chunks_size = 0
        for chunk in chunks:
            if current_chunks_size + (chunk[2] - chunk[1]) > chunk_size or len(
                    current_chunks) > max_nb_contig_per_chunk:
                current_chunks = []
                current_chunks_size = 0
                final_chunks.append(current_chunks)
            current_chunks.append(chunk)
            current_chunks_size += chunk[2] - chunk[1]
        return final_chunks

    def split_genome_chromosomes(self, with_unmapped=False):
        """
        Split the genome per chromosomes aggregating smaller chromosomes to similar length as the longest chromosome
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

    def split_genome_files(self, split_file_dir):
        """
        Write a bed file representing chunks of a reference genome, and return a dict where keys are the first chunk of
        each bed file and values are the file path:

        {
            ('chr1', 0, 15000000): 'path/to/dataset_region_chr1-0-15000000.bed'),
            ('chr2', 0, 20000000): 'path/to/dataset_region_chr2-0-20000000.bed')
        }
        """
        chunk_to_file = {}
        for c in self.genome_chunks:
            chunk_file = os.path.join(split_file_dir, self.dataset.name + '_region_%s-%s-%s.bed' % c[0])
            if not os.path.exists(chunk_file):
                with open(chunk_file, 'w') as open_file:
                    for chunk in c:
                        open_file.write('\t'.join([str(e) for e in chunk]) + '\n')
            chunk_to_file[c[0]] = chunk_file
        return chunk_to_file

    @staticmethod
    def chunks_from_fastq(indexed_fq_files):
        """
        Provide a list of chunks (start, end) for a pair of fastq file by parsing the index (.gbi) file.
        Each chunk contains 25M reads or less.
        """
        split_lines = 100000000
        grabix_index = indexed_fq_files[0] + '.gbi'
        assert os.path.isfile(grabix_index)
        with open(grabix_index) as open_file:
            next(open_file)  # discard first line of the file
            nb_lines = int(next(open_file))  # second line contains the number of lines in the file
        # Indices are 1 based end inclusive
        return split_in_chunks(nb_lines, split_lines, zero_based=False, end_inclusive=True)


class GATK4FilePath(Stage):
    """Provides hardcoded file path to subclasses"""

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
        return self.gatk4_basename + '_hard_filter_snps.vcf.gz'

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
        return self.dataset.genome_dict['data_files'].get('variation')

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
        vqsr_datasets = self.dataset.genome_dict['data_files'].get('vqsr')
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
    def vqsr_snps_r_script(self):
        return self.gatk4_basename + '_vqsr_snps.R'

    @property
    def vqsr_indels_output_recall(self):
        return self.gatk4_basename + '_vqsr_indels_recall.vcf.gz'

    @property
    def vqsr_indels_tranches(self):
        return self.gatk4_basename + '_vqsr_indels_tranches'

    @property
    def vqsr_indels_r_script(self):
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
        """Search existing fastq file for a specified run element."""
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
            os.path.join(
                self.indexed_fastq_dir, run_element.get(ELEMENT_RUN_ELEMENT_ID) + fastq_file[-len('_R1_001.fastq.gz'):]
            )
            for fastq_file in self._find_fastqs_for_run_element(run_element)
        ]

    @property
    def split_alignment_dir(self):
        return os.path.join(self.job_dir, 'split_alignment')

    def chunked_bam_file(self, run_element, chunk):
        return os.path.join(
            self.split_alignment_dir, run_element.get(ELEMENT_RUN_ELEMENT_ID) + '_name_sort_%s_%s.bam' % chunk
        )


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
            cpus=1,
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
                working_dir=self.exec_dir,
                cpus=1,
                mem=8
            ).join()
        return exit_status


class SplitBWA(SplitFastqStage):
    """Run bwa on chunks of fastq file provided by SplitFastqStage."""

    chunk_handler = Parameter()

    def bwa_command(self, fastq_pair, reference, expected_output_bam, read_group, chunk):
        tmp_file = expected_output_bam
        read_group_str = r'@RG\t%s' % r'\t'.join(['%s:%s' % (k, read_group[k]) for k in sorted(read_group)])

        command_bwa = '{bwa} mem -K 100000000 -Y -R \'{read_group}\' -M -t 2 {ref} ' \
                      '<({grabix} grab {fastq1} {chunk}) ' \
                      '<({grabix} grab {fastq2} {chunk})'
        command_bwa = command_bwa.format(bwa=toolset['bwa'], read_group=read_group_str, ref=reference,
                                         grabix=toolset['grabix'], fastq1=fastq_pair[0], fastq2=fastq_pair[1],
                                         chunk='%s %s' % (chunk[0], chunk[1]))
        alt_file = reference + '.alt'
        if os.path.isfile(alt_file):
            hla_out = os.path.splitext(expected_output_bam)[0] + '.hla'
            command_bwa += ' | {k8} {postalt} -p {hla_out} {alt_file}'.format(
                k8=toolset['k8'],
                postalt=toolset['postalt'],
                hla_out=hla_out,
                alt_file=alt_file
            )

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
                for chunk in self.chunk_handler.chunks_from_fastq(indexed_fq_files):
                    commands.append(self.bwa_command(
                        fastq_pair=indexed_fq_files,
                        reference=self.dataset.reference_genome,
                        expected_output_bam=self.chunked_bam_file(run_element, chunk),
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
    """Merge bam file generated in bwa mem and defined in SplitFastqStage."""

    chunk_handler = Parameter()

    def all_bam_chunks_file(self):
        all_chunk_bams = []
        for run_element in self.dataset.run_elements:
            if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                indexed_fq_files = self._indexed_fastq_for_run_element(run_element)
                for chunk in self.chunk_handler.chunks_from_fastq(indexed_fq_files):
                    all_chunk_bams.append(self.chunked_bam_file(run_element, chunk))
        bam_list_file = os.path.join(self.job_dir, self.dataset.name + '_all_bam_files.list')
        with open(bam_list_file, 'w') as open_file:
            open_file.write('\n'.join(all_chunk_bams))
        return bam_list_file

    def merge_command(self):
        cat_cmd = '{bamcat} level=0 tmpfile={cat_tmp} `cat {bam_list_file}`'.format(
            bamcat=toolset['bamcat'],
            cat_tmp=os.path.join(self.create_tmp_dir(), self.dataset.name),
            bam_list_file=self.all_bam_chunks_file()
        )
        bamsormadup_cmd = '{bamsormadup} threads=16 tmpfile={dup_tmp} indexfilename={merged_bam}.bai > ' \
                          '{merged_bam}'.format(
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
    """Generic class providing functions for running GATK commands"""
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

    def gatk_cmd(self, run_cls, output=None, memory=2, input=None, reference='default', spark_core=1, ext=None):
        return self._gatk_cmd(run_cls, output=output, memory=memory, input=input, reference=reference,
                              spark_core=spark_core, ext=ext)

    def gatk_picard_cmd(self, run_cls, output, memory=2, input=None, reference='default', ext=None):
        return self._gatk_cmd(run_cls, output=output, memory=memory, input=input, reference=reference,
                              ext=ext, picard_tool=True)


class PostAlignmentScatter(GATK4Stage):
    """Generic class providing ability to split the genome in chunks of approximately 20Mb"""
    @property
    def split_file_dir(self):
        d = os.path.join(self.job_dir, 'post_alignment_split')
        os.makedirs(d, exist_ok=True)
        return d

    def vcf_per_chunk(self, chunk):
        return os.path.join(self.split_file_dir, self.dataset.name + '_haplotype_caller_%s-%s-%s.vcf.gz' % chunk)


class SplitHaplotypeCaller(PostAlignmentScatter):
    """Run HaplotypeCaller on each chunk of genomes to create a VCF file"""

    bam_file = Parameter()
    chunk_handler = Parameter()

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
              for chunk, region_file in self.chunk_handler.split_genome_files(self.split_file_dir).items()],
            job_name='split_hap_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        return haplotype_status


class GatherVCF(PostAlignmentScatter):
    """Collate all vcf chunks into one"""

    chunk_handler = Parameter()

    def _run(self):
        vcf_list = os.path.join(self.split_file_dir, self.dataset.name + '_vcf.list')
        with open(vcf_list, 'w') as open_file:
            for chunks in self.chunk_handler.genome_chunks:
                open_file.write(self.vcf_per_chunk(chunks[0]) + '\n')

        concat_vcf_status = executor.execute(
            self.gatk_picard_cmd('GatherVcfs', self.genotyped_vcf, input=vcf_list, memory=6, reference=None),
            job_name='gather_geno_call',
            working_dir=self.exec_dir,
            cpus=1,
            mem=6
        ).join()

        if concat_vcf_status == 0:
            concat_vcf_status = common.tabix_vcf(self.exec_dir, self.genotyped_vcf)

        return concat_vcf_status


class SelectSNPs(GATK4Stage):
    """Extract the SNPs from a provided vcf file."""

    input_vcf = Parameter()

    def _run(self):
        select_var_command = self.gatk_cmd('SelectVariants', self.raw_snps_vcf)
        select_var_command += ' -V ' + self.input_vcf
        select_var_command += ' -select-type SNP '
        select_variants_status = executor.execute(
            select_var_command,
            job_name='snp_select',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()
        return select_variants_status


class SelectIndels(GATK4Stage):
    """Extract the Indels from a provided vcf file."""

    input_vcf = Parameter()

    def _run(self):
        select_var_command = self.gatk_cmd('SelectVariants', self.raw_indels_vcf)
        select_var_command += ' -V ' + self.input_vcf
        select_var_command += ' -select-type INDEL '
        select_variants_status = executor.execute(
            select_var_command,
            job_name='indel_select',
            working_dir=self.exec_dir,
            cpus=1,
            mem=16
        ).join()
        return select_variants_status


class SNPsFiltration(GATK4Stage):
    """Filter the SNPs using generic Hard filter."""

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
            cpus=1,
            mem=8
        ).join()
        return variant_filter_status


class IndelsFiltration(GATK4Stage):
    """Filter the Indels using generic Hard filter."""

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
            cpus=1,
            mem=16
        ).join()
        return variant_filter_status


class MergeVariants(GATK4Stage):
    """Merge multiple vcf files into one."""

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
            job_name='merge_vcf',
            working_dir=self.exec_dir,
            cpus=1,
            mem=8
        ).join()


class QCGATK4(Pipeline):
    toolset_type = 'gatk4_sample_processing'
    name = 'qc_gatk4'

    def __init__(self, dataset):
        super().__init__(dataset)
        self.chunk_handler = ChunkHandler(self.dataset)

    def build(self):
        """Build the QC pipeline."""
    
        # merge fastq to do contamination check
        merge_fastqs = self.stage(common.MergeFastqs)
        contam = self.stage(qc.FastqScreen, previous_stages=[merge_fastqs], fq_pattern=merge_fastqs.fq_pattern)
        blast = self.stage(qc.Blast, previous_stages=[merge_fastqs], fastq_file=merge_fastqs.fq_pattern.replace('?', '1'))
    
        # create fastq index then align via scatter-gather strategy
        fastq_index = self.stage(FastqIndex)
        split_bwa = self.stage(SplitBWA, previous_stages=[fastq_index])
        merge_bam_dup = self.stage(MergeBamAndDup, previous_stages=[split_bwa])
    
        # bam file QC
        samtools_stat = self.stage(common.SamtoolsStats, bam_file=merge_bam_dup.sorted_bam, previous_stages=[merge_bam_dup])
        samtools_depth = self.stage(qc.SamtoolsDepth, bam_file=merge_bam_dup.sorted_bam, previous_stages=[merge_bam_dup])
    
        # variants call via scatter-gather strategy
        haplotype_caller = self.stage(SplitHaplotypeCaller, bam_file=merge_bam_dup.sorted_bam, previous_stages=[merge_bam_dup])
        gather_vcf = self.stage(GatherVCF, previous_stages=[haplotype_caller])
    
        # variant filtering
        select_snps = self.stage(SelectSNPs, input_vcf=gather_vcf.genotyped_vcf, previous_stages=[gather_vcf])
        select_indels = self.stage(SelectIndels, input_vcf=gather_vcf.genotyped_vcf, previous_stages=[gather_vcf])
        filter_snps = self.stage(SNPsFiltration, previous_stages=[select_snps])
        filter_indels = self.stage(IndelsFiltration, previous_stages=[select_indels])
    
        # final variant merge
        merge_vcf = self.stage(MergeVariants,
                               vcf_files=[filter_snps.hard_filtered_snps_vcf, filter_indels.hard_filtered_indels_vcf],
                               output_vcf_file=filter_indels.hard_filtered_vcf,
                               previous_stages=[filter_snps, filter_indels])
    
        # variant file QC
        vcfstats = self.stage(qc.VCFStats, vcf_file=merge_vcf.hard_filtered_vcf, previous_stages=[merge_vcf])
    
        final_stages = [contam, blast, vcfstats, samtools_depth, samtools_stat]
    
        output = self.stage(common.SampleDataOutput, previous_stages=final_stages, output_fileset='gatk4_qc')
        review = self.stage(common.SampleReview, previous_stages=[output])
        cleanup = self.stage(common.Cleanup, previous_stages=[review])
    
        return cleanup
