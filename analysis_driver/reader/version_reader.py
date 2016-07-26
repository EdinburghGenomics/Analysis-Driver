import subprocess
import yaml
from analysis_driver.config import default as cfg


class VersionReader:
    def __init__(self, tool_name, command):
        self.tool_name = tool_name
        self.command = command

    def _get_stdout_from_command(self):
        p = subprocess.Popen(
            self.command.format(executable=cfg['tools'][self.tool_name], toolname=self.tool_name),
            stdout=subprocess.PIPE,
            shell=True
        )
        stdout, stderr = p.communicate()
        return stdout.strip().decode('utf-8')

    def _get_version_from_config(self):
        return cfg.query('versions', self.tool_name)

    def get_version(self):
        if self.command:
            return self._get_stdout_from_command()
        else:
            return self._get_version_from_config()


grep_version = '{executable} 2>&1 | grep "Version" | cut -d " " -f 2'  # e.g. 'Version: 0.1.2'
grep_toolname = '{executable} -v 2>&1 | grep "{toolname}" | cut -d " " -f 2 | head -n1'
executable_commands = {
    'bwa': grep_version,
    'samtools': grep_version,
    'bcl2fastq': grep_toolname,
    'fastqc': '{executable} -v 2>&1  | cut -d " " -f 2 ',
    'bcbio': '{executable} -v',
    'seqtk': grep_version,
    'samblaster': '{executable} -h 2>&1 | grep "Version" | cut -d " " -f 3',
    'sambamba': grep_toolname,
    'bamtools': grep_toolname,
    'bcftools': grep_toolname,
    'tabix': grep_version,
    'bgzip': grep_version,
    'fastqscreen': '{executable} -v 2>&1 | grep "fastq_screen" | cut -d " " -f 2 | head -n1',
    'sickle': '{executable} --version | grep "{toolname}" | cut -d " " -f 3 | head -n1',
    'verifybamid': '{executable} 2>&1 | grep "verifyBamID" | cut -d " " -f 2 | head -n1',
    'well_duplicate': None,
    'gatk': 'java -jar {executable} -h 2>&1 | grep "The Genome Analysis Toolkit (GATK)" | cut -d " " -f 6 | cut -d "," -f 1',
}


def get_versions():
    all_versions = {}
    for tool in cfg['tools']:
        if tool in executable_commands:
            all_versions[tool] = VersionReader(tool, executable_commands.get(tool)).get_version()
    return all_versions


def write_versions_to_yaml(yaml_file):
    with open(yaml_file, 'w') as o:
        o.write(yaml.safe_dump(get_versions(), default_flow_style=False))
