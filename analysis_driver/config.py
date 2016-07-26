import os.path
from egcg_core.config import Configuration, cfg


def _etc_config(config_file):
    return os.path.join(os.path.dirname(os.path.abspath(os.path.dirname(__file__))), 'etc', config_file)


def load_config():
    cfg.load_config_file(
        os.getenv('ANALYSISDRIVERCONFIG'),
        os.path.expanduser('~/.analysisdriver.yaml'),
        env_var='ANALYSISDRIVERENV'
    )

# backward compatibility
default = cfg

output_files_config = Configuration(_etc_config('output_files.yaml'))
sample_sheet_config = Configuration(_etc_config('sample_sheet_cfg.yaml'))
