from os.path import join
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.util import bash_commands
from analysis_driver.reader import RunInfo


class TestBashCommands(TestAnalysisDriver):
    def test_bcl2fastq(self):
        run_info = RunInfo(TestAnalysisDriver.assets_path)
        mask = run_info.reads.generate_mask(barcode_len=8)
        sample_sheet_csv = join(self.assets_path, 'SampleSheet_analysis_driver.csv')
        obs = bash_commands.bcl2fastq(self.assets_path, self.fastq_path, sample_sheet_csv, mask)
        exp = ('path/to/bcl2fastq_1.0.4 -l INFO --runfolder-dir %s --output-dir %s -r 8 -p 8 -w 8 '
               '--sample-sheet %s --use-bases-mask %s') % (self.assets_path, self.fastq_path, sample_sheet_csv, mask)
        assert obs == exp

    def test_fastqc(self):
        test_fastq = join(self.fastq_path, '10015AT', '10015ATA0001L05', 'this.fastq.gz')
        assert bash_commands.fastqc(test_fastq) == 'path/to/fastqc_v0.11.5 --nogroup -t 1 -q ' + test_fastq

    def test_bcbio(self):
        cmd = bash_commands.bcbio('run.yaml', self.assets_path)
        assert cmd == 'path/to/bcbio/bin/bcbio_nextgen.py run.yaml -n 10 --workdir ' + self.assets_path

    def test_prepare_samples(self):
        obs = bash_commands.bcbio_prepare_samples('a_job_dir', 'samples.csv')
        exp = 'path/to/bcbio/bin/bcbio_prepare_samples.py --out a_job_dir/merged --csv samples.csv'
        assert obs == exp

    def test_rsync_from_to(self):
        cmd = bash_commands.rsync_from_to('a_source', 'a_dest', exclude='an_exclude')
        assert cmd == 'rsync -rLD --exclude=an_exclude --update a_source a_dest'

    def test_bwa_mem_samblaster(self):
        cmd = bash_commands.bwa_mem_samblaster(['this.fastq', 'that.fastq'], 'ref.fasta', 'output/out.bam')
        expected_cmd = (
            'set -o pipefail; '
            'path/to/bwa_1.1 mem -M -t 16 ref.fasta this.fastq that.fastq | '
            'path/to/samblaster | '
            'path/to/samtools_1.3.1 view -b - | '
            'path/to/sambamba sort -m 5G --tmpdir output -t 16 -o output/out.bam /dev/stdin'
        )
        assert cmd == expected_cmd

    def test_bwa_mem_samblaster_read_groups(self):
        cmd = bash_commands.bwa_mem_samblaster(
            ['this.fastq', 'that.fastq'],
            'ref.fasta',
            'output/out.bam',
            read_group={'ID': '1', 'SM': 'user_sample_id', 'PL': 'illumina'}
        )
        expected_cmd = (
            'set -o pipefail; '
            'path/to/bwa_1.1 mem -M -t 16 -R \'@RG\\tID:1\\tPL:illumina\\tSM:user_sample_id\' '
            'ref.fasta this.fastq that.fastq | '
            'path/to/samblaster | '
            'path/to/samtools_1.3.1 view -b - | '
            'path/to/sambamba sort -m 5G --tmpdir output -t 16 -o output/out.bam /dev/stdin'
        )
        assert cmd == expected_cmd

    def test_bwa_mem_biobambam_read_groups(self):
        cmd = bash_commands.bwa_mem_biobambam(
            ['this.fastq', 'that.fastq'],
            'ref.fasta',
            'output/out.bam',
            read_group={'ID': '1', 'SM': 'user_sample_id', 'PL': 'illumina'}
        )
        expected_cmd = (
            'set -o pipefail; '
            'path/to/bwa_1.1 mem -M -t 16 -R \'@RG\\tID:1\\tPL:illumina\\tSM:user_sample_id\' '
            'ref.fasta this.fastq that.fastq | '
            'path/to/sortmapdup inputformat=sam SO=coordinate tmpfile=output/out.bam '
            'threads=16 indexfilename=output/out.bam.bai > output/out.bam'
        )
        assert cmd == expected_cmd

    def test_samtools_stats(self):
        expected = 'path/to/samtools_1.3.1 stats in.bam > out.txt'
        assert bash_commands.samtools_stats('in.bam', 'out.txt') == expected

    def test_samtools_depth_command(self):
        expected = 'path/to/samtools_1.3.1 depth -a -a -q 0 -Q 0 /path/to/bam_file | '\
                   'awk -F "	" \'{array[$1"	"$3]+=1} END{for (val in array){print val"	"array[val]}}\' | '\
                   'sort -T /path/to/job -k 1,1 -nk 2,2 > /path/to/depth_file'
        assert bash_commands.samtools_depth_command('/path/to/job', '/path/to/bam_file', '/path/to/depth_file') == expected

    def test_md5sum(self):
        assert bash_commands.md5sum('in.txt') == 'path/to/md5sum in.txt > in.txt.md5'

    def test_export_env_vars(self):
        cmds = bash_commands.export_env_vars()
        assert cmds == (
            'export PATH=path/to/bcbio/bin:path/to/jdk/bin:$PATH',
            'export LD_LIBRARY_PATH=path/to/bcbio/lib:$LD_LIBRARY_PATH',
            'export PERL5LIB=path/to/bcbio/lib/perl5:$PERL5LIB',
            'export JAVA_HOME=path/to/jdk',
            'export JAVA_BINDIR=path/to/jdk/bin',
            'export JAVA_ROOT=path/to/jdk',
            ''
        )

    def test_is_remote_path(self):
        assert not bash_commands.is_remote_path('a_file_path')
        assert bash_commands.is_remote_path('user@server:/home/user/a_file_path')

    def test_seqtk_fqchk(self):
        fastq_file = 'path/to/fastq_R1.fastq.gz'
        expected = 'path/to/seqtk fqchk -q 0 %s > %s.fqchk' % (fastq_file, fastq_file)
        assert bash_commands.seqtk_fqchk(fastq_file) == expected

    def test_fastq_filterer(self):
        assert bash_commands.fastq_filterer(('RE_R1_001.fastq.gz', 'RE_R2_001.fastq.gz'), 'RE_phix_read_name.txt') == (
            'run_filterer in_place RE_R1_001.fastq.gz RE_R2_001.fastq.gz RE_R1_001_filtered.fastq.gz '
            'RE_R2_001_filtered.fastq.gz RE_R1_001_filtered.fastq RE_R2_001_filtered.fastq RE_fastqfilterer.stats '
            'RE_phix_read_name.txt RE_R1_001.fastq_discarded RE_R2_001.fastq_discarded'
        )
        assert bash_commands.fastq_filterer(('RE_R1_001.fastq.gz', 'RE_R2_001.fastq.gz'), 'RE_phix_read_name.txt', trim_r1=149) == (
            'run_filterer keep_originals RE_R1_001.fastq.gz RE_R2_001.fastq.gz RE_R1_001_filtered.fastq.gz '
            'RE_R2_001_filtered.fastq.gz RE_R1_001_filtered.fastq RE_R2_001_filtered.fastq RE_fastqfilterer.stats '
            'RE_phix_read_name.txt RE_R1_001.fastq_discarded RE_R2_001.fastq_discarded --trim_r1 149'
        )

    def test_picard_gc_bias(self):
        cmd = bash_commands.picard_gc_bias(
            'directory/test.bam',
            'directory/metrics.txt',
            'directory/summary_metrics.txt',
            'directory/chart.pdf',
            'a_reference_genome.fasta'
        )

        assert cmd == ('path/to/java_8 -Djava.io.tmpdir=directory -XX:+UseSerialGC -Xmx8G -jar path/to/picard '
                       'CollectGcBiasMetrics INPUT=directory/test.bam OUTPUT=directory/metrics.txt '
                       'VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true CHART=directory/chart.pdf '
                       'R=a_reference_genome.fasta SUMMARY_OUTPUT=directory/summary_metrics.txt')

    def test_picard_mark_dup_command(self):
        cmd = bash_commands.picard_mark_dup_command(
            'directory/test.bam',
            'directory/test_out.bam',
            'directory/test_out.metrics'
        )
        exp = ('path/to/java_8 -Djava.io.tmpdir={tmpdir} -XX:+UseSerialGC -Xmx{mem}G -jar path/to/picard '
               'MarkDuplicates INPUT=directory/test.bam OUTPUT=directory/test_out.bam '
               'VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true METRICS_FILE=directory/test_out.metrics '
               'OPTICAL_DUPLICATE_PIXEL_DISTANCE=100')
        assert cmd == exp.format(tmpdir='directory', mem=10)

        cmd = bash_commands.picard_mark_dup_command(
            'directory/test.bam',
            'directory/test_out.bam',
            'directory/test_out.metrics',
            tmp_dir='directory2',
            memory=20
        )
        assert cmd == exp.format(tmpdir='directory2', mem=20)

    def test_picard_insert_size_command(self):
        cmd = bash_commands.picard_insert_size_command(
            'directory/test.bam',
            'directory/test_out.metrics',
            'directory/test_out.histogram'
        )
        exp = ('path/to/java_8 -Djava.io.tmpdir={tmpdir} -XX:+UseSerialGC -Xmx{mem}G -jar path/to/picard '
               'CollectInsertSizeMetrics INPUT=directory/test.bam OUTPUT=directory/test_out.metrics '
               'VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true HISTOGRAM_FILE=directory/test_out.histogram')
        assert cmd == exp.format(tmpdir='directory', mem=8)

        cmd = bash_commands.picard_insert_size_command(
            'directory/test.bam',
            'directory/test_out.metrics',
            'directory/test_out.histogram',
            tmp_dir='directory2',
            memory=10
        )
        assert cmd == exp.format(tmpdir='directory2', mem=10)

    def test_java_command(self):
        obs = bash_commands.java_command(1, 'a_tmp_dir', 'a_jar_file')
        assert obs == 'path/to/java_8 -Djava.io.tmpdir=a_tmp_dir -XX:+UseSerialGC -Xmx1G -jar a_jar_file '

    def test_tabix_command(self):
        assert bash_commands.tabix_vcf_command('file.vcf.gz') == 'path/to/tabix -f -p vcf file.vcf.gz'

    def test_bgzip_command(self):
        assert bash_commands.bgzip_command('file.vcf') == 'path/to/bgzip -f file.vcf'
