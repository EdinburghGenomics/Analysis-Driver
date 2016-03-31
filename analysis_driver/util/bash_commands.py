import os.path
from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg


app_logger = log_cfg.get_logger(__name__)


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


def seqtk_fqchk(fastq_file):
    cmd = cfg['tools']['seqtk'] + ' fqchk -q 0 %s > %s.fqchk' % (fastq_file, fastq_file)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def sickle_paired_end_in_place(fastq_file_pair):
    """
    Run sickle in paired end mode to do a very minimal trimming and filter reads shorter than 36 bases, i.e.
    remove the adapter dimers flagged by bcl2fastq.
    :param list fastq_file_pair: Two fastqs
    """
    if len(fastq_file_pair) != 2:
        raise AnalysisDriverError('sickle_paired_end only supports paired fastq files')

    f1, f2 = sorted(fastq_file_pair)
    name, ext = os.path.splitext(f1)
    of1 = name + '_sickle' + ext
    ofs = name + '_sickle_single' + ext
    lf = name + '_sickle.log'
    name, ext = os.path.splitext(f2)
    of2 = name + '_sickle' + ext
    cmds = []
    cmd = cfg['tools']['sickle'] + ' pe -f %s -r %s -o %s -p %s -s %s -q 5  -l 36  -x  -g -t sanger > %s'
    cmds.append(cmd % (f1, f2, of1, of2, ofs, lf))
    # replace the original files with the new files and remove the the single file to keep things clean
    cmds.append('EXIT_CODE=$?')
    cmds.append('(exit $EXIT_CODE) && mv %s %s' % (of1, f1))
    cmds.append('(exit $EXIT_CODE) && mv %s %s' % (of2, f2))
    cmds.append('(exit $EXIT_CODE) && rm %s' % ofs)
    cmds.append('(exit $EXIT_CODE)')
    for c in cmds:
        app_logger.debug('Writing: ' + c)
    return '\n'.join(cmds)


def bwa_mem_samblaster(fastq_pair, reference, expected_output_bam, thread=16):
    bwa_bin = cfg.query('tools', 'bwa')
    tmp_dir = os.path.dirname(expected_output_bam)
    command_bwa = '%s mem -M -t %s %s %s' % (bwa_bin, thread, reference, ' '.join(fastq_pair))
    command_samtools = '%s view -b -' % (cfg.query('tools', 'samtools'))
    command_sambamba = '%s sort -m 5G --tmpdir %s -t %s -o %s /dev/stdin' % (
        cfg.query('tools', 'sambamba'), tmp_dir, thread, expected_output_bam
    )
    cmd = ' | '.join([command_bwa, cfg.query('tools', 'samblaster'), command_samtools, command_sambamba])
    app_logger.debug('Writing: ' + cmd)
    return cmd


def bamtools_stats(bam_file, output_file):
    bamtools_bin = cfg.query('tools', 'bamtools')
    cmd = '%s stats -in %s -insert > %s' % (bamtools_bin, bam_file, output_file)
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


def rsync_from_to(source, dest, exclude=None):
    command = 'rsync -ruLD '
    if exclude:
        command += '--exclude=%s ' % exclude
    if is_remote_path(source) or is_remote_path(dest):
        command += '-e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -c arcfour" '

    command += '%s %s' % (source, dest)
    return command
