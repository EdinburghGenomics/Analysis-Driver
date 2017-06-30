from os.path import join
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.util import bash_commands
from analysis_driver.reader import RunInfo

assets = TestAnalysisDriver.assets_path
fastq_path = TestAnalysisDriver.fastq_path


def test_bcl2fastq():
    run_info = RunInfo(TestAnalysisDriver.assets_path)
    mask = run_info.reads.generate_mask(barcode_len=8)
    sample_sheet_csv = join(assets, 'SampleSheet_analysis_driver.csv')
    obs = bash_commands.bcl2fastq(assets, fastq_path, sample_sheet_csv, mask)
    exp = ('path/to/bcl2fastq_1.0.4 -l INFO --runfolder-dir %s --output-dir %s -r 8 -d 8 -p 8 -w 8 '
           '--sample-sheet %s --use-bases-mask %s') % (assets, fastq_path, sample_sheet_csv, mask)
    assert obs == exp


def test_fastqc():
    test_fastq = join(fastq_path, '10015AT', '10015ATA0001L05', 'this.fastq.gz')
    expected = '--nogroup -t 1 -q ' + test_fastq
    assert bash_commands.fastqc(test_fastq).endswith(expected)


def test_bcbio():
    cmd = bash_commands.bcbio('run.yaml', assets)
    assert cmd == 'path/to/bcbio/bin/bcbio_nextgen.py run.yaml -n 10 --workdir ' + assets


def test_prepare_samples():
    obs = bash_commands.bcbio_prepare_samples('a_job_dir', 'samples.csv')
    exp = 'path/to/bcbio/bin/bcbio_prepare_samples.py --out a_job_dir/merged --csv samples.csv'
    assert obs == exp


def test_rsync_from_to():
    cmd = bash_commands.rsync_from_to('a_source', 'a_dest', exclude='an_exclude')
    assert cmd == 'rsync -rLD --exclude=an_exclude --update a_source a_dest'


def test_bwa_mem_samblaster():
    cmd = bash_commands.bwa_mem_samblaster(['this.fastq', 'that.fastq'], 'ref.fasta', 'output/out.bam')
    expected_cmd = (
        'set -o pipefail; '
        'path/to/bwa_1.1 mem -M -t 16 ref.fasta this.fastq that.fastq | '
        'path/to/samblaster | '
        'path/to/samtools view -b - | '
        'path/to/sambamba sort -m 5G --tmpdir output -t 16 -o output/out.bam /dev/stdin'
    )
    assert cmd == expected_cmd


def test_bwa_mem_samblaster_read_groups():
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
        'path/to/samtools view -b - | '
        'path/to/sambamba sort -m 5G --tmpdir output -t 16 -o output/out.bam /dev/stdin'
    )
    assert cmd == expected_cmd


def test_bwa_mem_biobambam_read_groups():
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
        'path/to/bamsormadup inputformat=sam SO=coordinate tmpfile=output/out.bam '
        'threads=16 indexfilename=output/out.bam.bai > output/out.bam'
    )
    assert cmd == expected_cmd


def test_samtools_stats():
    expected = 'path/to/samtools stats in.bam > out.txt'
    assert bash_commands.samtools_stats('in.bam', 'out.txt') == expected


def test_md5sum():
    assert bash_commands.md5sum('in.txt') == 'md5sum in.txt > in.txt.md5'


def test_export_env_vars():
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


def test_is_remote_path():
    assert not bash_commands.is_remote_path('a_file_path')
    assert bash_commands.is_remote_path('user@server:/home/user/a_file_path')


def test_seqtk_fqchk():
    fastq_file = 'path/to/fastq_R1.fastq.gz'
    expected = 'path/to/seqtk fqchk -q 0 %s > %s.fqchk' % (fastq_file, fastq_file)
    assert bash_commands.seqtk_fqchk(fastq_file) == expected


def test_fastq_filterer_and_pigz_in_place():
    exp = (
        'run_filterer RE_R1_001.fastq.gz RE_R2_001.fastq.gz RE_R1_001_filtered.fastq.gz '
        'RE_R2_001_filtered.fastq.gz RE_R1_001_filtered.fastq RE_R2_001_filtered.fastq '
        '--stats_file RE_fastqfilterer.stats'
    )
    assert bash_commands.fastq_filterer_and_pigz_in_place(('RE_R1_001.fastq.gz', 'RE_R2_001.fastq.gz')) == exp
