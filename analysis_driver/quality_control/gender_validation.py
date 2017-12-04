import os
from egcg_core import executor, util
from analysis_driver.segmentation import Parameter, Stage


class GenderValidation(Stage):
    vcf_file = Parameter()

    def _run(self):
        """Detect gender of the sample based on the %het on the X chromosome."""
        name, ext = os.path.splitext(util.find_file(self.vcf_file))
        if ext == '.gz':
            file_opener = 'zcat'
            name, gz = os.path.splitext(name)
        else:
            file_opener = 'cat'

        gender_call_file = name + '.sex'

        command = util.str_join(
            '%s %s' % (file_opener, self.vcf_file),
            "grep -P '^chrX|^X'",
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
