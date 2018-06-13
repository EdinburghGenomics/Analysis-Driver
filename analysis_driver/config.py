import os.path
from egcg_core.config import Configuration, cfg


def etc_config(config_file):
    return os.path.join(os.path.dirname(os.path.abspath(os.path.dirname(__file__))), 'etc', config_file)


def load_config():
    cfg.load_config_file(
        os.getenv('ANALYSISDRIVERCONFIG'),
        os.path.expanduser('~/.analysisdriver.yaml'),
        env_var='ANALYSISDRIVERENV'
    )


class OutputFileConfiguration(Configuration):
    def __init__(self):
        super().__init__(etc_config('output_files.yaml'))
        self.pipeline_set = False

    def set_pipeline_type(self, pipeline_type):
        if not self.pipeline_set:
            self.content = self.content[pipeline_type]
            self.pipeline_set = True

    def job_dir_file(self, outfile_id):
        outfile_record = self.content.get(outfile_id)
        if outfile_record:
            location = os.path.join(*outfile_record.get('location', ['']))
            return os.path.join(location, outfile_record['basename'])

    def output_dir_file(self, outfile_id):
        outfile_record = self.content.get(outfile_id)
        if outfile_record:
            return outfile_record.get('new_name') or outfile_record['basename']


default = cfg  # backward compatibility
sample_sheet_config = Configuration(etc_config('sample_sheet_cfg.yaml'))
tool_versioning_cfg = Configuration(etc_config('tool_versioning.yaml'))
output_file_config = OutputFileConfiguration()