__author__ = 'mwham'
import sys
sys.path.append('..')
import analysis_driver
sys.exit(analysis_driver.driver.main(sys.argv[1:]))
