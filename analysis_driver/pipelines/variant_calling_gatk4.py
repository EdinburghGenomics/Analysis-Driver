import os
from os.path import join

from egcg_core import executor, util
from egcg_core.config import cfg
from egcg_core.constants import ELEMENT_NB_READS_CLEANED, ELEMENT_RUN_ELEMENT_ID, ELEMENT_PROJECT_ID, ELEMENT_RUN_NAME, \
    ELEMENT_LANE

from analysis_driver import segmentation
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.tool_versioning import toolset
from analysis_driver.util.bash_commands import java_command

toolset_type = 'gatk4_sample_processing'
name = 'variant_calling_gatk4'


class SplitFastqStage(segmentation.Stage):

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
            os.path.join(self.indexed_fastq_dir, run_element.get(ELEMENT_RUN_ELEMENT_ID) + fastq_file[-len('_R1.fastq.gz'):])
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

    def _run(self):
        """get all fastq files, uncompress and recompress, then index with grabix."""
        commands = []
        for run_element in self.dataset.run_elements:
            if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                for ifastqs, ofastqs in zip(
                        self._find_fastqs_for_run_element(run_element),
                        self._indexed_fastq_for_run_element(run_element)
                ):
                    cmd_template = 'gunzip -c {ipath} | {pbgzip} -n 16  -c /dev/stdin > {opath}'
                    commands.append(cmd_template.format(ipath=ifastqs[0], pbgzip=toolset['pbgzip'], opath=ofastqs[0]))
                    commands.append(cmd_template.format(ipath=ifastqs[1], pbgzip=toolset['pbgzip'], opath=ofastqs[1]))
        exit_status1 = executor.execute(
            commands,
            job_name='compress_fastq',
            working_dir=self.job_dir
        ).join()
        exit_status2 = 0
        if exit_status1 != 0:
            commands = []
            for run_element in self.dataset.run_elements:
                if int(run_element.get(ELEMENT_NB_READS_CLEANED, 0)) > 0:
                    for fastq_file in self._indexed_fastq_for_run_element(run_element):
                        commands.append(grabix=toolset['grabix'] + ' index ' + fastq_file)
            exit_status2 = executor.execute(
                commands,
                job_name='index_fastq',
                working_dir=self.job_dir
            ).join()
        return exit_status1 + exit_status2


class SplitBWA(SplitFastqStage):

    def bwa_command(self, fastq_pair, reference, expected_output_bam, read_group, chunk):
        tmp_file = expected_output_bam
        read_group_str = r'@RG\t%s' % r'\t'.join(['%s:%s' % (k, read_group[k]) for k in sorted(read_group)])

        command_bwa = '{bwa} -R \'{read_group}\' mem -M  {ref} ' \
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
        """get all fastq files, and align them in chunks."""

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
            cpu=1,
            memory=2
        ).join()
        return exit_status


class MergeBam(SplitFastqStage):

    def all_chunk_bam_list_file(self):
        all_chunk_bams = []
        for run_element in self.dataset.run_elements:
            indexed_fq_files = self._indexed_fastq_for_run_element(run_element)
            for chunk in self.chunks_from_fastq(indexed_fq_files):
                all_chunk_bams.append(self.chuncked_bam_file(run_element, chunk))

    def merge_command(self):
        cat_cmd = '{bamcat} level=0 tmpfile={cat_tmp} `cat {bam_list_file}`'.format(
            bamcat=toolset['bamcat'],
            cat_tmp='',
            bam_list_file=self.all_chunk_bam_list_file()
        )
        bamsormadup_cmd = '{bamsormadup} threads=16 tmpfile={dup_tmp} indexfilename={merged_bam}.bai > {merged_bam}'.format(
            bamsormadup=toolset['bamsormadup'],
            dup_tmp='',
            merged_bam=''
        )
        return 'set -o pipefail; ' + ' | '.join([cat_cmd, bamsormadup_cmd])

    def _run(self):
        exit_status = executor.execute(
            self._merge_command(),
            job_name='merge_dup_bam',
            working_dir=self.job_dir,
            cpu=12,
            memory=24
        ).join()
        return exit_status


class GATK4Stage(segmentation.Stage):
    @property
    def gatk_run_dir(self):
        d = os.path.join(self.job_dir, 'gatk4')
        os.makedirs(d, exist_ok=True)
        return d

    def gatk_cmd(self, run_cls, output, memory=2, input_bam=None, cpu=16, ext=None):
        base_cmd = toolset['gatk_bin'] + \
                   ' --java-options "-Djava.io.tmpdir={tmp_dir} -XX:+UseSerialGC -Xmx{memory}G" ' + \
                   '{run_cls} -R {ref} -o {output} -l INFO ' \
                   '--spark-master local[{cpu}] ' \
                   '--conf spark.driver.host=localhost ' \
                   '--conf spark.network.timeout=800 ' \
                   '--conf spark.local.dir={tmp_dir} '
        if input_bam:
            base_cmd += ' -I {input_bam}'
        if ext:
            base_cmd += ext

        return base_cmd.format(
            tmp_dir=self.gatk_run_dir,
            cpu=cpu,
            memory=memory,
            ref=self.dataset.reference_genome,
            run_cls=run_cls,
            input_bam=input_bam,
            output=output
        )

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
    def indel_realigned_bam(self):
        return self.basename + '_indel_realigned.bam'

    @property
    def sample_gvcf(self):
        return self.basename + '.g.vcf'

    @property
    def genotyped_vcf(self):
        return self.basename + '.vcf'

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


class BaseRecal(GATK4Stage):
    def _run(self):
        '''
        java
        -Dsamjdk.use_async_io_read_samtools=false
        -Dsamjdk.use_async_io_write_samtools=true
        -Dsamjdk.use_async_io_write_tribble=false
        -Dsamjdk.compression_level=2
        -Xms500m -Xmx45864m
        -Djava.io.tmpdir=tmpodVBCD
        -jar gatk-package-4.0.3.0-local.jar
        BaseRecalibratorSpark -I 104619_2010_A-sort.bam
        --spark-master local[16]
        --output tmpodVBCD/104619_2010_A-sort-recal.grp
        --reference hg38.2bit
        --conf spark.driver.host=localhost
        --conf spark.network.timeout=800
        --conf spark.local.dir=tmpodVBCD
        --known-sites dbsnp-150.vcf.gz
        -L cleaned-104619_2010_A-sort-callable_sample.bed
        --interval-set-rule INTERSECTION
        '''
        return executor.execute(
            self.gatk_cmd(
                'HaplotypeCallerSpark', self.output_grp, input_bam=self.sorted_bam,
                xmx=48, spark_core=16, ext=' --knownSites ' + self.dbsnp
            ),
            job_name='gatk_base_recal',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()


class PrintReads(GATK4Stage):
    def _run(self):
        return executor.execute(
            self.gatk_cmd(
                'PrintReads', self.recal_bam, input_bam=self.sorted_bam,
                xmx=48, nt=1, ext=' -BQSR ' + self.output_grp
            ),
            job_name='gatk_print_reads',
            working_dir=self.gatk_run_dir,
            cpus=16,
            mem=64
        ).join()



def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)
    fastq_index = stage(FastqIndex)
    split_bwa = stage(SplitBWA, previous_stages=[fastq_index])
    merge_bam = stage(MergeBam, previous_stages=[split_bwa])
    return merge_bam
