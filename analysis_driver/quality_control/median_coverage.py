import os
from threading import Thread
from analysis_driver.app_logging import AppLogger
from analysis_driver.config import default as cfg

class GatkDepthofCoverage(AppLogger, Thread):
    def __init__(self, bam_file):
        self.bam_file = bam_file
        Thread.__init__(self)

    def _get_gatk_depthofcoverage_command(self):
        reference = cfg['bwa'] + 'genomes/Hsapiens/' + cfg['genome'] + '/seq/' + cfg['genome'] + '.fa'






    def run(self):
        try:
            self._get_median_coverage()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.median_coverage_expected_outfiles