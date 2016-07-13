from os.path import join, dirname, abspath
from egcg_core.config import cfg
from analysis_driver.config import default

cfg.load_config_file(default.config_file)
__version__ = join(dirname(dirname(abspath(__file__))), 'version.txt').strip().lstrip('v')
