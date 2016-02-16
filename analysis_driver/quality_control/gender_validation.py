import argparse
import os
from threading import Thread
import sys
from analysis_driver.app_logging import AppLogger
from analysis_driver import executor, util
from analysis_driver.config import default as cfg


class GenderValidation(AppLogger, Thread):
    """
    This class will perform the Gender validation steps. It subclasses Thread, allowing it to run in the
    background.
    """
    def __init__(self, sample_id, vcf_file):
        self.vcf_file = vcf_file
        self.sample_id = sample_id
        self.exception = None
        self.return_value = None
        Thread.__init__(self)

    def _gender_call(self):
        """
        Detect gender of the sample based on the %het on the X chromosome.
        :rtype: list
        :return list of file containing the results of the validation.
        """

        name, ext = os.path.splitext(self.vcf_file)
        if ext == '.gz':
            file_opener = 'zcat'
            name, dummy = os.path.splitext(name)
        else:
            file_opener = 'cat'

        gender_call_file = name + '.sex'

        command = util.str_join(
            '%s %s' % (file_opener, self.vcf_file),
            "grep '^chrX'",
            "awk '{split($10,a,\":\"); count[a[1]]++; total++} END{for (g in count){print g\" \"count[g]/total}}'",
            "grep '0/1'",
            "awk '{if ($2>.35){gender=\"FEMALE\"}else{if ($2<.15){gender=\"MALE\"}else{gender=\"UNKNOWN\"}} print gender}'",
            separator=' | '
        ) + ' > ' + gender_call_file
        self.info(command)

        return executor.execute(
            [command],
            job_name='sex_detection',
            run_id=self.sample_id,
            walltime=6,
            cpus=1,
            mem=2,
            log_command=False
        ).join()

    def run(self):
        try:
            self.return_value = self._gender_call()
        except Exception as e:
            self.exception = e

    def join(self, timeout=None):
        super().join(timeout=timeout)
        if self.exception:
            raise self.exception
        return self.return_value


def main():
    args = _parse_args()
    os.makedirs(os.path.join(cfg['jobs_dir'], args.sample_id), exist_ok=True)
    s = GenderValidation(args.sample_id, args.vcf_file)
    s.start()
    return s.join()


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-v', '--vcf_file', dest="vcf_file", type=str, help='the vcf file used to detect the gender')
    p.add_argument('-s', '--sample_id', dest="sample_id", type=str, help='the sample id to be used as job directory')

    return p.parse_args()

if __name__ == "__main__":
    sys.exit(main())
