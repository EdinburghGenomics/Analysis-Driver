__author__ = 'mwham'
import sys
import os.path

sys.path.append(os.path.dirname(os.path.dirname(__file__)))

import analysis_driver
sys.exit(analysis_driver.driver.main(sys.argv[1:]))
