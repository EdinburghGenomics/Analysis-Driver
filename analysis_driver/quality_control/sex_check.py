import os
from egcg_core import executor, util
from analysis_driver.segmentation import Parameter, Stage

_sex_aliases = {'female': ['f', 'female', 'girl', 'woman'], 'male': ['m', 'male', 'boy', 'man']}


def sex_alias(sex):
    for key in _sex_aliases:
        if str(sex).lower() in _sex_aliases[key]:
            return key
    return 'unknown'


class SexCheck(Stage):
    vcf_file = Parameter()

    def _run(self):
        """Detect sex of the sample based on the %het on the X chromosome."""
        name, ext = os.path.splitext(util.find_file(self.vcf_file))
        if ext == '.gz':
            file_opener = 'zcat'
            name, gz = os.path.splitext(name)
        else:
            file_opener = 'cat'

        sex_check_file = name + '.sex'

        command = util.str_join(
            '%s %s' % (file_opener, self.vcf_file),
            "grep -P '^chrX|^X'",
            "awk '{split($10,a,\":\"); count[a[1]]++; total++} END{for (g in count){print g\" \"count[g]/total}}'",
            "grep '0/1'",
            "awk '{if ($2>.35){sex=\"FEMALE\"}else{if ($2<.15){sex=\"MALE\"}else{sex=\"UNKNOWN\"}} print sex, $2}'",
            separator=' | '
        ) + ' > ' + sex_check_file
        self.info(command)

        return executor.execute(
            command,
            job_name='sex_check',
            working_dir=self.job_dir,
            walltime=6,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()
