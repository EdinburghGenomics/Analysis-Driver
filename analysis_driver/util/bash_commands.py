import os.path
from analysis_driver.config import default as cfg
from egcg_core.app_logging import logging_default as log_cfg
from analysis_driver.exceptions import AnalysisDriverError


app_logger = log_cfg.get_logger('bash_commands')


def bcl2fastq(input_dir, fastq_path, sample_sheet=None, mask=None):
    """
    Build a bcl2fastq command for an Illumina HiSeqX dataset.
    :param input_dir: Path to the input dir containing the bcl files
    :param fastq_path: Path to the dir in which to send generated .fastqs
    :param str sample_sheet: Path to a bcl2fast2 sample sheet to use
    :param str mask: A mask to use, as generated by reader.RunInfo
    """
    cmd = '%s -l INFO --runfolder-dir %s --output-dir %s -r 8 -d 8 -p 8 -w 8' % (
        cfg['tools']['bcl2fastq'], input_dir, fastq_path
    )
    if sample_sheet:
        cmd += ' --sample-sheet ' + sample_sheet
    if mask:
        cmd += ' --use-bases-mask ' + mask
    app_logger.debug('Writing: ' + cmd)
    return cmd


def fastqc(fastq, threads=1):
    cmd = cfg['tools']['fastqc'] + ' --nogroup -t %s -q %s' % (threads, fastq)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def gzip_test(f):
    cmd = 'gzip -t ' + f
    app_logger.debug('Writing: ' + cmd)
    return cmd


def seqtk_fqchk(fastq_file):
    cmd = cfg['tools']['seqtk'] + ' fqchk -q 0 %s > %s.fqchk' % (fastq_file, fastq_file)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def fq_filt_prelim_cmd():
    cmd = (
        'function run_filterer {{',  # double up { and } to escape them for str.format()
        'i1=$1', 'i2=$2', 'o1=$3', 'o2=$4', 'fifo_1=$5', 'fifo_2=$6', 'log_file=$7',

        'mkfifo $fifo_1', 'mkfifo $fifo_2',
        '{ff} --i1 $i1 --i2 $i2 --o1 $fifo_1 --o2 $fifo_2 --threshold {threshold} $* > $log_file 2>&1 &',
        'fq_filt_pid=$!',
        '{pigz} -c -p {pzt} $fifo_1 > $o1 &',
        'pigz_r1_pid=$!',
        '{pigz} -c -p {pzt} $fifo_2 > $o2 &',
        'pigz_r2_pid=$!',

        'exit_status=0',
        'wait $fq_filt_pid',
        'exit_status=$[$exit_status + $?]',
        'wait $pigz_r1_pid',
        'exit_status=$[$exit_status + $?]',
        'wait $pigz_r2_pid',
        'exit_status=$[$exit_status + $?]',
        'rm $fifo_1 $fifo_2',
        'if [ $exit_status == 0 ]; then mv $o1 $i1; mv $o2 $i2; fi',
        '(exit $exit_status)',
        '}}'
    )
    return '\n'.join(cmd).format(
        ff=cfg['tools']['fastq-filterer'],
        threshold=cfg.query('fastq_filterer', 'min_length', ret_default='36'),
        pigz=cfg.query('tools', 'pigz', ret_default='pigz'),
        pzt=10  # two pigz processes run, so pigz threads here will be doubled
    )


def fastq_filterer_and_pigz_in_place(fastq_file_pair, tiles_to_filter=None, trim_r1=None, trim_r2=None):
    """
    :param tuple[str,str] fastq_file_pair: Paired-end fastqs to filter
    :param list tiles_to_filter: Tile IDs for reads to remove regardless of length
    :param trim_r1: Maximum length to trim R1 to
    :param trim_r2: As trim_r1, but for R2
    Run fastq filterer on a pair of fastqs, removing pairs where one is shorter than 36 bases
    """
    if len(fastq_file_pair) != 2:
        raise AnalysisDriverError('fastq-filterer only supports paired fastq files')

    i1, i2 = sorted(fastq_file_pair)

    base_1 = i1.replace('.fastq.gz', '')
    fifo_1 = base_1 + '_filtered.fastq'
    o1 = fifo_1 + '.gz'

    base_2 = i2.replace('.fastq.gz', '')
    fifo_2 = base_2 + '_filtered.fastq'
    o2 = fifo_2 + '.gz'

    stats_file = base_1.replace('_R1_001', '') + '_fastqfilterer.stats'
    log_file = base_1.replace('_R1_001', '') + '.log'

    fastq_filterer_cmd = 'run_filterer {0} {1} {2} {3} {4} {5} {6} --stats_file {7}'.format(
        i1, i2, o1, o2, fifo_1, fifo_2, log_file, stats_file
    )

    if tiles_to_filter:
        fastq_filterer_cmd += ' --remove_tiles %s' % (','.join([str(t) for t in tiles_to_filter]))
    if trim_r1:
        fastq_filterer_cmd += ' --trim_r1 %s' % trim_r1
    if trim_r2:
        fastq_filterer_cmd += ' --trim_r2 %s' % trim_r2

    return fastq_filterer_cmd


def bwa_mem_samblaster(fastq_pair, reference, expected_output_bam, read_group=None, thread=16):
    tmp_dir = os.path.dirname(expected_output_bam)
    command_bwa = '%s mem -M -t %s' % (cfg['tools']['bwa'], thread)

    if read_group:
        read_group_str = '@RG\\t%s' % '\\t'.join(['%s:%s' % (k, read_group[k]) for k in sorted(read_group)])
        command_bwa += ' -R \'%s\'' % read_group_str

    command_bwa += ' %s %s' % (reference, ' '.join(fastq_pair))
    command_samtools = cfg['tools']['samtools'] + ' view -b -'
    command_sambamba = '%s sort -m 5G --tmpdir %s -t %s -o %s /dev/stdin' % (
        cfg['tools']['sambamba'], tmp_dir, thread, expected_output_bam
    )
    cmd = 'set -o pipefail; ' + ' | '.join([command_bwa, cfg['tools']['samblaster'], command_samtools, command_sambamba])
    app_logger.debug('Writing: ' + cmd)
    return cmd


def bwa_mem_biobambam(fastq_pair, reference, expected_output_bam, read_group=None, thread=16):
    tmp_file = expected_output_bam
    index = expected_output_bam + '.bai'
    command_bwa = '%s mem -M -t %s' % (cfg['tools']['bwa'], thread)

    if read_group:
        read_group_str = '@RG\\t%s' % '\\t'.join(['%s:%s' % (k, read_group[k]) for k in sorted(read_group)])
        command_bwa += ' -R \'%s\'' % read_group_str

    command_bwa += ' %s %s' % (reference, ' '.join(fastq_pair))
    command_bambam = '%s inputformat=sam SO=coordinate tmpfile=%s threads=%s indexfilename=%s > %s'
    command_bambam = command_bambam % (cfg['tools']['biobambam_sortmapdup'], tmp_file, thread, index, expected_output_bam)

    cmd = 'set -o pipefail; ' + ' | '.join([command_bwa, command_bambam])
    app_logger.debug('Writing: ' + cmd)
    return cmd


def samtools_stats(bam_file, output_file):
    samtools_bin = cfg.query('tools', 'samtools')
    cmd = '%s stats %s > %s' % (samtools_bin, bam_file, output_file)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def md5sum(input_file):
    cmd = cfg.query('tools', 'md5sum', ret_default='md5sum') + ' %s > %s.md5' % (input_file, input_file)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def export_env_vars():
    """Write export statements for environment variables required by BCBio"""
    app_logger.debug('Writing Java paths')
    return (
        _export('PATH', os.path.join(cfg['tools']['bcbio'], 'bin'), prepend=True),
        _export('PATH', os.path.join(cfg['tools']['jdk'], 'bin'), prepend=True),
        _export('LD_LIBRARY_PATH', os.path.join(cfg['tools']['bcbio'], 'lib'), prepend=True),
        _export('PERL5LIB', os.path.join(cfg['tools']['bcbio'], 'lib', 'perl5'), prepend=True),
        _export('JAVA_HOME', cfg['tools']['jdk']),
        _export('JAVA_BINDIR', os.path.join(cfg['tools']['jdk'], 'bin')),
        _export('JAVA_ROOT', cfg['tools']['jdk']),
        ''
    )


def bcbio(run_yaml, workdir, threads=10):
    cmd = '%s %s -n %s --workdir %s' % (
        os.path.join(cfg['tools']['bcbio'], 'bin', 'bcbio_nextgen.py'), run_yaml, threads, workdir
    )
    app_logger.debug('Writing: ' + cmd)
    return cmd


def _export(env_var, value, prepend=False):
    statement = 'export %s=%s' % (env_var, value)
    if prepend:
        statement += ':$' + env_var
    return statement


def is_remote_path(fp):
    return (':' in fp) and ('@' in fp)


def rsync_from_to(source, dest, exclude=None, size_only=False):
    command = 'rsync -rLD '
    if exclude:
        command += '--exclude=%s ' % exclude
    if size_only:
        command += '--size-only '
    else:
        command += '--update '
    if is_remote_path(source) or is_remote_path(dest):
        command += '-e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -c arcfour" '

    command += '%s %s' % (source, dest)
    return command


def java_command(memory, tmp_dir, jar):
    return 'java -Djava.io.tmpdir={tmp_dir} -XX:+UseSerialGC -Xmx{memory}G -jar {jar}'.format(
        memory=memory,
        tmp_dir=tmp_dir,
        jar=jar
    )
