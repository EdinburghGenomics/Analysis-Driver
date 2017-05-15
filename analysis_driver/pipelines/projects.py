import os
from egcg_core.util import find_file
from analysis_driver.config import default as cfg, OutputFileConfiguration
from analysis_driver.quality_control import Relatedness
from analysis_driver.exceptions import PipelineError
from analysis_driver.transfer_data import output_project_data
from analysis_driver import segmentation


def build_pipeline(dataset):
    project_source = os.path.join(cfg.query('sample', 'delivery_source'), dataset.name)
    gvcf_files = []
    for sample in dataset.samples_processed:
        gvcf_file = find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
        if gvcf_file:
            gvcf_files.append(gvcf_file)
    if len(gvcf_files) < 2:
        raise PipelineError('Incorrect number of gVCF files: require at least two')

    relatedness = Relatedness(
        dataset=dataset,
        gvcf_files=gvcf_files,
        reference=cfg['references'][dataset.species]['fasta'],
        project_id=dataset.name
    )
    output = Output(dataset=dataset, previous_stages=[relatedness])
    return output


class Output(segmentation.Stage):
    def _run(self):

        dir_with_output_files = os.path.join(self.job_dir, 'relatedness_outfiles')
        os.makedirs(dir_with_output_files, exist_ok=True)
        output_files_cfg = OutputFileConfiguration('project_process')

        for symlink_file in output_files_cfg.content.values():
            source = os.path.join(
                self.job_dir,
                os.path.join(*symlink_file['location']),
                symlink_file['basename'].format(project_id=self.dataset.name)
            )

            symlink_path = os.path.join(dir_with_output_files, symlink_file['basename'].format(project_id=self.dataset.name))
            if os.path.isfile(source):
                if os.path.islink(symlink_path):
                    os.unlink(symlink_path)
                os.symlink(source, symlink_path)
            else:
                raise PipelineError('Could not find the file ' + source + ', unable to create link')

        return output_project_data(dir_with_output_files, self.dataset.name)
