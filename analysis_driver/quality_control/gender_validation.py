import argparse
import os
import sys
from luigi import Parameter
from egcg_core import executor, util
from analysis_driver.segmentation import Stage


class GenderValidation(Stage):
    vcf_file = Parameter()

    def _run(self):
        """Detect gender of the sample based on the %het on the X chromosome."""
        name, ext = os.path.splitext(self.vcf_file)
        if ext == '.gz':
            file_opener = 'zcat'
            name, gz = os.path.splitext(name)
        else:
            file_opener = 'cat'

        gender_call_file = name + '.sex'

        command = util.str_join(
            '%s %s' % (file_opener, self.vcf_file),
            "grep '^chrX'",
            "awk '{split($10,a,\":\"); count[a[1]]++; total++} END{for (g in count){print g\" \"count[g]/total}}'",
            "grep '0/1'",
            "awk '{if ($2>.35){gender=\"FEMALE\"}else{if ($2<.15){gender=\"MALE\"}else{gender=\"UNKNOWN\"}} print gender, $2}'",
            separator=' | '
        ) + ' > ' + gender_call_file
        self.info(command)

        return executor.execute(
            command,
            job_name='sex_detection',
            working_dir=self.job_dir,
            walltime=6,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


def main():
    from analysis_driver.config import default as cfg
    from analysis_driver.dataset_scanner import SampleScanner
    args = _parse_args()
    os.makedirs(args.working_dir, exist_ok=True)
    dataset = SampleScanner(cfg).get_dataset(args.sample_id)
    s = GenderValidation(dataset=dataset, vcf_file=args.vcf_file)
    s.start()
    return s.join()


def _parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--sample_id', type=str, help='sample ID for creating a Sample dataset object')
    p.add_argument('-v', '--vcf_file', dest='vcf_file', type=str, help='the vcf file used to detect the gender')
    p.add_argument('-s', '--working_dir', dest='working_dir', type=str, help='the working dir for execution')
    return p.parse_args()

if __name__ == '__main__':
    sys.exit(main())
