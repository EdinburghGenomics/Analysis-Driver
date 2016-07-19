import subprocess
import yaml

from analysis_driver.exceptions import AnalysisDriverError
from analysis_driver.config import default as cfg


class Version_reader():
    tool_name = None
    command = None
    def __init__(self, tool_name, command):
        self.tool_name = tool_name
        self.command = command

    def _get_stdout_from_command(self):
        if self.command:
            tool = cfg.query('tools', self.tool_name)
            p = subprocess.Popen(tool + self.command, stdout=subprocess.PIPE)
            stdout = p.stdout.read()
            p.wait()
            p.stdout.close()
            return stdout

    def get_version_from_config(self):
        return cfg.query('versions', self.tool_name)

    def get_version(self):
        version = self._get_stdout_from_command()
        if not version:
            version = self.get_version_from_config()
        if version:
            return version
        return None

command1 = ' 2>&1 | grep "Version" | cut -d " " -f 2'

all_commands = {
    'bwa': command1,
    'samtools': command1,
    'bcl2fastq': ' -v 2>&1 | grep "bcl2fastq" | cut -d " " -f 2 | head -n1 ',
    'fastqc': ' -v 2>&1  | cut -d " " -f 2 ',
    'bcbio': ' -v',
    'seqtk': command1,
    'samblaster': ' -h 2>&1 | grep "Version" | cut -d " " -f 3',
    'sambamba': ' -v 2>&1 | grep "sambamba" | cut -d " " -f 2 | head -n1',
    'bamtools': ' -v 2>&1 | grep "bamtools" | cut -d " " -f 2 | head -n1',
    'gatk': ' -h 2>&1 | grep "The Genome Analysis Toolkit (GATK)" | cut -d " " -f 6 | cut -d "," -f 1',
    'bcftools': ' -v 2>&1 | grep "bcftools" | cut -d " " -f 2 | head -n1',
    'tabix': command1,
    'fastqscreen':' -v 2>&1 | grep "fastq_screen" | cut -d " " -f 2 | head -n1',
    'sickle': ' --version | grep "sickle" | cut -d " " -f 3 | head -n1',
    'verifybamid': ' 2>&1 | grep "verifyBamID" | cut -d " " -f 2 | head -n1',
    'well_duplicate': None
}


def get_versions():
    all_versions = {}
    for tool in cfg['tools']:
        v = Version_reader(tool, all_commands.get(tool)).get_version()
        if v :
            all_versions[tool] = v
    return all_versions

def write_versions_to_yaml(yaml_file):
    with open(yaml_file, 'w') as o:
        o.write(yaml.safe_dump(get_versions(), default_flow_style=False))