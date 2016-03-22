__author__ = 'mwham'
import os.path
from tests.test_analysisdriver import TestAnalysisDriver
from analysis_driver.util import bash_commands
from analysis_driver.reader import SampleSheet, RunInfo
from analysis_driver.config import default as cfg

helper = TestAnalysisDriver()
sample_sheet_csv = os.path.join(helper.assets_path, 'SampleSheet_analysis_driver.csv')
sample_sheet = SampleSheet(helper.sample_sheet_path)
run_info = RunInfo(helper.assets_path)


def test_bcl2fastq():
    mask = sample_sheet.generate_mask(run_info.mask)
    helper.compare_lists(
        observed=bash_commands.bcl2fastq(helper.assets_path, helper.fastq_path, sample_sheet_csv, mask),
        expected=(
            cfg['tools']['bcl2fastq'] +
            ' -l INFO'
            ' --runfolder-dir ' + helper.assets_path +
            ' --output-dir ' + helper.fastq_path +
            ' -r 8 -d 8 -p 8 -w 8' +
            ' --sample-sheet ' + sample_sheet_csv +
            ' --use-bases-mask ' + mask
        )
    )


def test_fastqc():
    test_fastq = os.path.join(helper.fastq_path, '10015AT', '10015ATA0001L05', 'this.fastq.gz')
    expected = '--nogroup -t 1 -q ' + test_fastq
    assert bash_commands.fastqc(test_fastq).endswith(expected)


def test_bcbio():
    cmd = bash_commands.bcbio('run.yaml', helper.assets_path)
    assert cmd == 'path/to/bcbio/bin/bcbio_nextgen.py run.yaml -n 10 --workdir ' + helper.assets_path


def test_rsync_from_to():
    cmd = bash_commands.rsync_from_to('a_source', 'a_dest', exclude='an_exclude')
    assert cmd == 'rsync -ruLD --exclude=an_exclude a_source a_dest'


def test_bwa_mem_samblaster():
    cmd = bash_commands.bwa_mem_samblaster(['this.fastq', 'that.fastq'], 'ref.fasta', 'output/out.bam')
    expected_cmd = (
        'path/to/bwa mem -M -t 16 ref.fasta this.fastq that.fastq | '
        'path/to/samblaster | '
        'path/to/samtools view -b - | '
        'path/to/sambamba sort -m 5G --tmpdir output -t 16 -o output/out.bam /dev/stdin'
    )
    print(cmd)
    print(expected_cmd)
    assert cmd == expected_cmd


def test_bamtools_stats():
    expected = 'path/to/bamtools stats -in in.bam -insert > out.txt'
    assert bash_commands.bamtools_stats('in.bam', 'out.txt') == expected


def test_md5sum():
    assert bash_commands.md5sum('in.txt') == 'md5sum in.txt > in.txt.md5'


def test_export_env_vars():
    cmds = bash_commands.export_env_vars()
    assert cmds == (
        'export PATH=path/to/bcbio/bin:$PATH',
        'export LD_LIBRARY_PATH=path/to/bcbio/lib:$LD_LIBRARY_PATH',
        'export PERL5LIB=path/to/bcbio/lib/perl5:$PERL5LIB',
        'export JAVA_HOME=path/to/jdk',
        'export JAVA_BINDIR=path/to/jdk/bin',
        'export JAVA_ROOT=path/to/jdk',
        ''
    )


def test_export():
    assert bash_commands._export('THIS', 'that') == 'export THIS=that'
    assert bash_commands._export('THIS', 'that', prepend=True) == 'export THIS=that:$THIS'


def test_is_remote_path():
    assert not bash_commands.is_remote_path('a_file_path')
    assert bash_commands.is_remote_path('user@server:/home/user/a_file_path')


def test_seqtk_fqchk():
    fastq_file = 'path/to/fastq_R1.fastq.gz'
    expected = '%s fqchk -q 0 %s > %s.fqchk' % (cfg['tools']['seqtk'], fastq_file, fastq_file)
    assert bash_commands.seqtk_fqchk(fastq_file) == expected
