__author__ = 'mwham'
import sys
import os


main_cfg = {
    
    'str': (
        ('clarity', 'baseuri'),
        ('clarity', 'username'),
        ('clarity', 'password'),
        'genome',
        'job_execution',
        'job_queue',
        'qsub',
        'bcl2fastq',
        'fastqc',
        'bcbio',
        'jdk'
    ),
    'dict': {
        'rest_api': ('url',),
        'logging': ('format', 'datefmt'),
        ('notification', 'email_notification'): ('mailhost', 'port', 'reporter_email', 'recipient_emails'),
        'genotype-validation': (
            'reference',
            'genotypes_repository',
            'bwa',
            'gatk',
            'samblaster',
            'samtools',
            'sambamba'
        ),
        'clarity': ('baseuri', 'username', 'password')
    }
}


def validate_main_config(cfg):

    for q in main_cfg['str']:
        _check(cfg, 'str', q)

    for q, keys in main_cfg.get('dict', {}).items():
        _check(cfg, 'dict', q, dict_keys=keys)

    invalid_file_paths = _validate_file_paths(cfg.content)
    if invalid_file_paths:
        warn('invalid file paths: ' + str(invalid_file_paths))

    for name, h in cfg.query('logging', 'handlers', ret_default={}).items():
        if 'stream' in h or 'filename' in h:
            pass
        else:
            warn('logging handler %s has %s.' % (name, str(h.keys())))
            print('should have either \'stream\' or \'filename\' specified')


def _check(cfg, param_type, query, dict_keys=None):
    if type(query) is str:
        query = (query,)

    check = None
    param = cfg.query(*query, ret_default=[])

    if param_type == 'file_path':
        check = os.path.isfile(param)
    elif param_type == 'dir':
        check = os.path.isdir(param)
    elif param_type == 'str':
        check = type(param) is str
    elif param_type == 'int':
        check = type(param) is int
    elif param_type == 'dict':
        check = all([k in param for k in dict_keys])
    elif param_type == 'list':
        check = type(param) is list
    else:
        warn('unknown type specified for param %s: \'%s\'' % (param, param_type))

    if check is False:
        warn('param %s (from %s) is not a valid %s' % (param, '/'.join(query), param_type))
    return check


def _validate_file_paths(content=None):
    """
    Recursively search through the values of self.content and if the value is an absolute file path,
    assert that it exists.
    :param content: a dict, list or str (i.e. potential file path) to validate
    """
    invalid_file_paths = []
    if type(content) is dict:
        vals = list(content.values())
        for v in vals:
            invalid_file_paths.extend(_validate_file_paths(v))
    elif type(content) is list:
        for v in content:
            invalid_file_paths.extend(_validate_file_paths(v))
    elif type(content) is str:
        if content.startswith('/') and not os.path.exists(content):
            invalid_file_paths.append(content)

    return invalid_file_paths


def warn(msg):
    print('warning: ' + msg)


if __name__ == '__main__':
    try:
        os.putenv('ANALYSISDRIVERCONFIG', sys.argv[1])
    except IndexError:
        pass

    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from analysis_driver import config

    validate_main_config(config.default)
