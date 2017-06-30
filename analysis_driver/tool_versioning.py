import subprocess
from egcg_core.app_logging import AppLogger
from analysis_driver.config import cfg, tool_versioning_cfg
from analysis_driver.exceptions import AnalysisDriverError


class Toolset(AppLogger):
    tools = {}
    tool_versions = {}
    toolset_version = None
    toolset_valid = True
    latest_version = max(tool_versioning_cfg['toolsets'].keys())

    def select_toolset_version(self, version):
        self.toolset_version = version
        self.tools = {}
        self.tool_versions = {}

        cfg_entry = tool_versioning_cfg['toolsets'][self.toolset_version]

        for k in cfg_entry:
            self.add_tool(k)

        if not self.toolset_valid:
            raise AnalysisDriverError('Invalid toolset for version %s' % self.toolset_version)
        else:
            self.info('Selected toolsett version %s', self.toolset_version)

    def add_tool(self, toolname):
        version = self.resolve(toolname, 'version')
        version_cmd = self.resolve(toolname, 'version_cmd')

        self.tools[toolname] = cfg['tools'][toolname][version]
        self.tool_versions[toolname] = version
        self.check_version(toolname, version, version_cmd)

    def check_version(self, toolname, version, version_cmd):
        if version_cmd in tool_versioning_cfg['cmd_aliases']:
            version_cmd = tool_versioning_cfg['cmd_aliases'][version_cmd]

        observed_version = self._get_stdout(
            version_cmd.format(executable=self.tools[toolname], toolname=toolname)
        )
        if observed_version != str(version):
            self.error('%s tool version: expected %s, got %s', toolname, version, observed_version)
            self.toolset_valid = False

    def resolve(self, toolname, key, source_toolset=None):
        if source_toolset is None:
            source_toolset = self.toolset_version
        elif source_toolset < 0:
            raise AnalysisDriverError('Could not resolve %s for %s', key, toolname)

        config = tool_versioning_cfg['toolsets'][source_toolset][toolname] or {}

        if key in config:
            self.debug('Resolved tool %s with %s %s from toolset %s', toolname, key, config[key], source_toolset)
            return config[key]
        else:
            return self.resolve(toolname, key, source_toolset - 1)

    @staticmethod
    def _get_stdout(cmd):
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
        return stdout.strip().decode('utf-8')

    def __getitem__(self, item):
        return self.tools[item]

toolset = Toolset()
