import os
from egcg_core.util import find_file
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.config import default as cfg, OutputFileConfiguration
from analysis_driver.quality_control import Relatedness, Peddy, Genotype_gVCFs
from analysis_driver.exceptions import PipelineError
from analysis_driver.transfer_data import output_project_data
from analysis_driver import segmentation


def build_pipeline(dataset):
    project_id = dataset.name
    sample_ids = [sample['sample_id'] for sample in dataset.samples_processed]

    project_source = os.path.join(cfg.query('sample', 'delivery_source'), project_id)
    gvcf_files = []
    for sample in dataset.samples_processed:
        gvcf_file = find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
        if not gvcf_file:
            raise PipelineError('Unable to find gVCF file for sample %s' % (sample))
        gvcf_files.append(gvcf_file)
    if len(gvcf_files) < 2:
        raise PipelineError('Incorrect number of gVCF files: require at least two')

    reference = cfg['references'][dataset.species]['fasta']
    genotype_gvcfs = Genotype_gVCFs(dataset=dataset, gVCFs=gvcf_files, reference=reference)
    relatedness = Relatedness(dataset=dataset, previous_stages=[genotype_gvcfs])
    peddy = Peddy(dataset=dataset, ids=sample_ids, previous_stages=[genotype_gvcfs])
    output = Output(dataset=dataset, previous_stages=[relatedness, peddy])
    cleanup = Cleanup(previous_stages=[output])
    return cleanup


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
