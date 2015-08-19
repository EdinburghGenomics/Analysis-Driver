__author__ = 'mwham'
import sys
import os.path

if __name__ == '__main__':

    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    import analysis_driver.process_trigger

    sys.exit(analysis_driver.process_trigger.main())
