import os
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.pipelines import common
from analysis_driver import segmentation


class DataOutput(segmentation.Stage):
    def _run(self):
        dir_with_linked_files = common.link_results_files(self.dataset.name, self.job_dir, 'non_human_qc')
        write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))
        return common.output_data(self.dataset, self.job_dir, self.dataset.name, dir_with_linked_files)


class Cleanup(segmentation.Stage):
    def _run(self):
        return common.cleanup(self.dataset.name)


def build_pipeline(dataset):
    bam_file_production = common.build_bam_file_production(dataset)
    data_output = DataOutput(dataset=dataset, previous_stages=bam_file_production)
    cleanup = Cleanup(dataset=dataset, previous_stages=[data_output])
    return cleanup
