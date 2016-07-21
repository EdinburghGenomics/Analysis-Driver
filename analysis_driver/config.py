import os.path
from egcg_core.config import Configuration, EnvConfiguration


def _etc_config(config_file):
    return os.path.join(os.path.dirname(os.path.abspath(os.path.dirname(__file__))), 'etc', config_file)


# singletons for access by other modules
default = EnvConfiguration(
    os.getenv('ANALYSISDRIVERCONFIG'),
    os.path.expanduser('~/.analysisdriver.yaml'),
    _etc_config('example_analysisdriver.yaml')
)
output_files_config = Configuration(_etc_config('output_files.yaml'))
sample_sheet_config = Configuration(_etc_config('sample_sheet_cfg.yaml'))
