__author__ = 'mwham'
import os.path
from analysis_driver.app_logging import get_logger
from analysis_driver.config import default as cfg


app_logger = get_logger('command_writer')


def bcl2fastq(mask, input_dir, fastq_path):
    """
    :param str mask: A mask to use, as generated by reader.RunInfo
    :param input_dir: Path to the input dir containing the bcl files
    :param fastq_path: Path to the dir in which to send generated .fastqs
    :rtype: str
    """
    cmd = ' '.join(['%s -l INFO --runfolder-dir %s',
                    '--output-dir %s --sample-sheet %s',
                    '--use-bases-mask %s',
                    '-r 8 -d 8 -p 8 -w 8']) % (
        cfg['bcl2fastq'],
        input_dir,
        fastq_path,
        os.path.join(input_dir, 'SampleSheet_analysis_driver.csv'),
        mask
    )

    app_logger.debug('Writing command: ' + cmd)
    return cmd


def fastqc(fastq, threads=4):
    """
    :param str fastq: An input fastq file
    :rtype: str
    """
    cmd = cfg['fastqc'] + ' --nogroup -t %s -q %s' % (threads, fastq)
    app_logger.debug('Writing: ' + cmd)
    return cmd


def bcbio_java_paths():
    """
    Write export statements for Java environment variables required by the execution of BCBio pipelines using
    GATK.
    :rtype: list[str]
    """
    cmds = []
    app_logger.debug('Writing Java paths')

    cmds.append('export JAVA_HOME=' + cfg['gatk']['java_home'])
    cmds.append('export JAVA_BINDIR=' + cfg['gatk']['java_bindir'])
    cmds.append('export JAVA_ROOT=' + cfg['gatk']['java_root'] + '\n')

    return cmds


def bcbio(run_yaml, workdir, cores=16):
    """
    :param run_yaml: The pipeline config file to be run
    :param workdir: The desired working directory
    :param cores: The number of cores for BCBio to use
    :rtype: str
    """
    cmd = ('%s %s -n %s --workdir %s' % (cfg['bcbio_nextgen'], run_yaml, cores, workdir))
    app_logger.debug('Writing: ' + cmd)
    return cmd
