import yaml
import subprocess
from os.path import isfile
from shutil import copyfile
from multiprocessing import Lock
from egcg_core.app_logging import AppLogger
from analysis_driver.config import cfg, tool_versioning_cfg
from analysis_driver.exceptions import AnalysisDriverError


class Toolset(AppLogger):
    def __init__(self):
        self.tools = {}
        self.tool_versions = {}
        self.versioning_cfg = tool_versioning_cfg
        self.version = None
        self.type = None
        self.versions_file = None
        self.lock = Lock()

    def latest_version(self, toolset_type):
        return max(self.versioning_cfg['toolsets'][toolset_type])

    def configure(self, toolset_type, toolset_version, versions_file):
        self.type = toolset_type
        self.version = toolset_version
        self.versions_file = versions_file

        self.tools = {}
        self.tool_versions = {}
        self.info('Selected %s toolset version %s', self.type, self.version)

    @property
    def unversioned_tools(self):
        return [k for k in self.tools if k not in self.tool_versions]

    def add_tool(self, toolname):
        if toolname in self.versioning_cfg['toolsets'][self.type][self.version]:
            self.tools[toolname] = self.add_versioned_tool(toolname)
        else:
            executable = cfg['tools'].get(toolname)
            if not executable:
                raise AnalysisDriverError('Referenced tool %s not present in config' % toolname)

            if isinstance(executable, list):
                executable = sorted(executable)[-1]

            self.tools[toolname] = executable

    def add_versioned_tool(self, toolname):
        """
        For a given tool, resolve its version and version_cmd from self.versioning_cfg and find the correct executable
        in cfg['tools'].
        :param str toolname:
        """

        version = self.resolve(toolname, 'version')
        version_cmd = self.resolve(toolname, 'version_cmd')
        dependency = self.resolve(toolname, 'dependency')
        if dependency:
            self.add_tool(dependency)

        executables = cfg['tools'][toolname]
        if type(executables) is str:
            executables = [executables]

        for e in executables:
            if self.check_version(toolname, e, version, version_cmd, self.tools.get(dependency)):
                self.tool_versions[toolname] = version
                self.update_versions_file(toolname, version)
                self.debug('Resolved %s version %s: %s', toolname, version, e)
                return e

        raise AnalysisDriverError('Could not find version %s for %s' % (version, toolname))

    def check_version(self, toolname, executable, exp_version, version_cmd, dependency=None):
        """
        Query an executable for its version and compare with an expected value.
        :param str toolname:
        :param str executable:
        :param str exp_version:
        :param str version_cmd:
        :param str dependency:
        """
        if version_cmd in self.versioning_cfg['cmd_aliases']:
            version_cmd = self.versioning_cfg['cmd_aliases'][version_cmd]

        obs_version = self._get_stdout(
            version_cmd.format(executable=executable, toolname=toolname, dependency=dependency)
        )
        return obs_version == str(exp_version)

    def resolve(self, toolname, val, toolset_version=None):
        """
        Extract a value from self.versioning_cfg for a given tool, inheriting from previous toolset versions if needed.
        :param str toolname:
        :param str val: The value to extract for the tool, e.g. version or version_cmd
        :param int toolset_version:
        """
        if toolset_version is None:
            toolset_version = self.version  # place to start scanning backwards to 0
        elif toolset_version < 0:  # finished scanning backwards and found nothing
            self.debug('%s not resolved for %s', val, toolname)
            return None

        config = self.versioning_cfg['toolsets'][self.type][toolset_version].get(toolname) or {}

        if val in config:
            return config[val]
        else:
            return self.resolve(toolname, val, toolset_version - 1)

    def write_to_yaml(self, file_path):
        copyfile(self.versions_file, file_path)

    def update_versions_file(self, toolname, version):
        with self.lock:
            if isfile(self.versions_file):
                with open(self.versions_file, 'r') as f:
                    tool_versions = yaml.safe_load(f)
            else:
                tool_versions = {}

            tool_versions[toolname] = version

            with open(self.versions_file, 'w') as f:
                f.write(yaml.safe_dump(tool_versions, default_flow_style=False))

    @staticmethod
    def _get_stdout(cmd):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        return stdout.strip().decode('utf-8')

    def __getitem__(self, item):
        if item not in self.tools:
            self.add_tool(item)

        return self.tools[item]


toolset = Toolset()
