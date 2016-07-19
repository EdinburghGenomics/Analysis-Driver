from os.path import join, dirname, abspath
__version__ = join(dirname(dirname(abspath(__file__))), 'version.txt').strip().lstrip('v')
