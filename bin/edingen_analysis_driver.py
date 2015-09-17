__author__ = 'mwham'
import sys
import os.path

if __name__ == '__main__':

    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import analysis_driver.client

    sys.exit(analysis_driver.client.main())
