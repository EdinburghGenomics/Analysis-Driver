import argparse
import os
from threading import Thread
import sys
from analysis_driver.app_logging import AppLogger
from analysis_driver import executor, util


class GenderValidation(AppLogger, Thread):
    """
    This class will perform the Gender validation steps. It subclasses Thread, allowing it to run in the
    background.
    """
    def __init__(self, working_dir, vcf_file):
        self.vcf_file = vcf_file
        self.working_dir = working_dir
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
            working_dir=self.working_dir,
            walltime=6,
            cpus=1,
            mem=2,
            log_commands=False
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
    os.makedirs(args.working_dir, exist_ok=True)
    s = GenderValidation(args.working_dir, args.vcf_file)
    s.start()
    return s.join()


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-v', '--vcf_file', dest="vcf_file", type=str, help='the vcf file used to detect the gender')
    p.add_argument('-s', '--working_dir', dest="working_dir", type=str, help='the working dir for execution')
    return p.parse_args()

if __name__ == "__main__":
    sys.exit(main())
