import subprocess
from egcg_core.app_logging import AppLogger
from analysis_driver.config import cfg, tool_versioning_cfg
from analysis_driver.exceptions import AnalysisDriverError


class Toolset(AppLogger):
    tools = {}
    tool_versions = {}
    toolset_version = None
    latest_version = max(tool_versioning_cfg['toolsets'].keys())

    def select_toolset(self, version):
        self.toolset_version = version
        self.tools = {}
        self.tool_versions = {}

        for k in cfg['tools']:
            if k in tool_versioning_cfg['toolsets'][self.toolset_version]:
                self.tools[k] = self.add_versioned_tool(k)
            else:
                self.tools[k] = cfg['tools'][k]

        self.info('Selected toolset version %s', self.toolset_version)

    def add_versioned_tool(self, toolname):
        version = self.resolve(toolname, 'version')
        version_cmd = self.resolve(toolname, 'version_cmd')

        config = cfg['tools'][toolname]
        if type(config) is str:
            config = [config]

        for c in config:
            if self.check_version(toolname, c, version, version_cmd):
                self.tool_versions[toolname] = version
                return c

        raise AnalysisDriverError('Could not find version %s for %s' % (version, toolname))

    def check_version(self, toolname, executable, exp_version, version_cmd):
        if version_cmd in tool_versioning_cfg['cmd_aliases']:
            version_cmd = tool_versioning_cfg['cmd_aliases'][version_cmd]

        obs_version = self._get_stdout(version_cmd.format(executable=executable, toolname=toolname))
        return obs_version == str(exp_version)

    def resolve(self, toolname, key, toolset_version=None):
        if toolset_version is None:
            toolset_version = self.toolset_version
        elif toolset_version < 0:
            raise AnalysisDriverError('Could not resolve %s for %s', key, toolname)

        config = tool_versioning_cfg['toolsets'][toolset_version][toolname] or {}

        if key in config:
            self.debug(
                'Resolved tool %s with %s %s from toolset %s', toolname, key, config[key], toolset_version
            )
            return config[key]
        else:
            return self.resolve(toolname, key, toolset_version - 1)

    @staticmethod
    def _get_stdout(cmd):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        return stdout.strip().decode('utf-8')

    def __getitem__(self, item):
        return self.tools[item]

toolset = Toolset()
