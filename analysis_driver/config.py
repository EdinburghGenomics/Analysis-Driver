import os.path
from egcg_core.config import Configuration, cfg


def _etc_config(config_file):
    return os.path.join(os.path.dirname(os.path.abspath(os.path.dirname(__file__))), 'etc', config_file)

default_search_path = [
    os.getenv('PROJECTMANAGEMENTCONFIG'),
    os.path.expanduser('~/.project_management.yaml')
]

def load_config(*config_files):
    if not config_files:
        config_files = default_search_path
    cfg.load_config_file(*config_files, env_var=os.getenv('ANALYSISDRIVERENV'))

#still have a default analysis driver config backward compatibility
default = cfg
output_files_config = Configuration(_etc_config('output_files.yaml'))
sample_sheet_config = Configuration(_etc_config('sample_sheet_cfg.yaml'))
