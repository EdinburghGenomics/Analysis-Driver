import yaml
import subprocess
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

    def latest_version(self, toolset_type):
        return max(self.versioning_cfg['toolsets'][toolset_type])

    def select_type(self, toolset_type):
        self.type = toolset_type

    def select_version(self, version):
        if self.type is None:
            raise AnalysisDriverError('Tried to select a toolset version with no type set')

        self.version = version
        self.tools = {}
        self.tool_versions = {}

        for k in cfg['tools']:
            self.add_tool(k)

        self.info('Selected %s toolset version %s', self.type, self.version)

    @property
    def unversioned_tools(self):
        return [k for k in self.tools if k not in self.tool_versions]

    def add_tool(self, toolname):
        if toolname in self.versioning_cfg['toolsets'][self.type][self.version]:
            self.tools[toolname] = self.add_versioned_tool(toolname)
        else:
            config = cfg['tools'][toolname]
            if isinstance(config, list):
                config = sorted(config)[-1]

            self.tools[toolname] = config

    def add_versioned_tool(self, toolname):
        """
        For a given tool, resolve its version and version_cmd from self.versioning_cfg and find the correct executable
        in cfg['tools'].
        :param str toolname:
        """
        if toolname in self.tools:
            self.debug('Tried to add %s - already added')
            return self.tools[toolname]

        version = self.resolve(toolname, 'version')
        version_cmd = self.resolve(toolname, 'version_cmd')
        dependency = self.resolve(toolname, 'dependency')
        if dependency:
            self.add_tool(dependency)

        config = cfg['tools'][toolname]
        if type(config) is str:
            config = [config]

        for c in config:
            if self.check_version(toolname, c, version, version_cmd, self.tools.get(dependency)):
                self.tool_versions[toolname] = version
                return c

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
            self.warning('Could not resolve %s for %s', val, toolname)
            return None

        config = self.versioning_cfg['toolsets'][self.type][toolset_version][toolname] or {}

        if val in config:
            self.debug(
                'Resolved tool %s with %s=%s from toolset %s', toolname, val, config[val], toolset_version
            )
            return config[val]
        else:
            return self.resolve(toolname, val, toolset_version - 1)

    def write_to_yaml(self, file_path):
        with open(file_path, 'w') as f:
            f.write(yaml.safe_dump(self.tool_versions, default_flow_style=False))

    @staticmethod
    def _get_stdout(cmd):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        return stdout.strip().decode('utf-8')

    def __getitem__(self, item):
        return self.tools[item]


toolset = Toolset()
