__author__ = 'mwham'
import os.path
from analysis_driver.app_logging import logging_default as log_cfg
from analysis_driver.config import default as cfg


app_logger = log_cfg.get_logger('command_writer')


def bcl2fastq(input_dir, fastq_path, sample_sheet=None, mask=None):
    """
    :param str mask: A mask to use, as generated by reader.RunInfo
    :param input_dir: Path to the input dir containing the bcl files
    :param fastq_path: Path to the dir in which to send generated .fastqs
    :rtype: str
    """
    cmd = '%s -l INFO --runfolder-dir %s --output-dir %s -r 8 -d 8 -p 8 -w 8' % (
        cfg['tools']['bcl2fastq'], input_dir, fastq_path
    )
    if sample_sheet:
        cmd += ' --sample-sheet ' + sample_sheet
    if mask:
        cmd += ' --use-bases-mask ' + mask
    app_logger.debug('Writing command: ' + cmd)
    return cmd


def fastqc(fastq, threads=1):
    cmd = cfg['tools']['fastqc'] + ' --nogroup -t %s -q %s' % (threads, fastq)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def seqtk_fqchk(fastq_file):
    cmd = cfg['tools']['seqtk'] + ' fqchk -q 0 %s > %s.fqchk' % (fastq_file, fastq_file)
    app_logger.debug('Writing: ' + cmd)
    return cmd


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
    """
    :param run_yaml: The generated yaml config file to be run for the pipeline
    """
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
    """rsync command that will transfer the file to the desired destination"""
    command = 'rsync -rLD '
    if exclude:
        command += '--exclude=%s ' % exclude
    if is_remote_path(source) or is_remote_path(dest):
        command += '-e "ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -c arcfour" '

    command += '%s %s' % (source, dest)
    return command
