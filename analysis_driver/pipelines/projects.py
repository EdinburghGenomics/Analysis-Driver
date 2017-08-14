import os
from egcg_core import executor
from egcg_core.util import find_file
from analysis_driver import segmentation
from analysis_driver.util import bash_commands
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.config import default as cfg, OutputFileConfiguration
from analysis_driver.quality_control import Relatedness, Peddy, GenotypeGVCFs, ParseRelatedness
from analysis_driver.exceptions import PipelineError
from analysis_driver.transfer_data import output_data_and_archive, create_output_links

toolset_type = 'project_processing'
name = 'project'


def delivery_source():
    return cfg.query('sample', 'delivery_source')

def build_pipeline(dataset):
    sample_ids = [sample['sample_id'] for sample in dataset.samples_processed]
    project_source = os.path.join(delivery_source(), dataset.name)
    gvcf_files = []
    for sample in dataset.samples_processed:
        gvcf_file = find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
        if not gvcf_file:
            raise PipelineError('Unable to find gVCF file for sample %s in %s' % (sample['sample_id'], project_source))
        gvcf_files.append(gvcf_file)
    if len(gvcf_files) < 2:
        raise PipelineError('Incorrect number of gVCF files: require at least two')

    reference = cfg['references'][dataset.species]['fasta']
    genotype_gvcfs = GenotypeGVCFs(dataset=dataset, gVCFs=gvcf_files, reference=reference)
    relatedness = Relatedness(dataset=dataset, previous_stages=[genotype_gvcfs])
    peddy = Peddy(dataset=dataset, ids=sample_ids, previous_stages=[genotype_gvcfs])
    parse = ParseRelatedness(dataset=dataset, ids=sample_ids, parse_method='parse_both', previous_stages=[relatedness, peddy])
    md5 = MD5Sum(dataset=dataset, previous_stages=[parse])
    output = Output(dataset=dataset, previous_stages=[md5])
    cleanup = Cleanup(dataset=dataset, previous_stages=[output])
    return cleanup


class MD5Sum(segmentation.Stage):
    def _run(self):
        dir_with_linked_files = os.path.join(self.job_dir, 'relatedness_outfiles')
        return executor.execute(
            *[bash_commands.md5sum(os.path.join(dir_with_linked_files, f)) for f in os.listdir(dir_with_linked_files)],
            job_name='md5sum',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()


class Output(segmentation.Stage):
    def _run(self):
        dir_with_linked_files = os.path.join(self.job_dir, 'relatedness_outfiles')
        os.makedirs(dir_with_linked_files, exist_ok=True)
        output_cfg = OutputFileConfiguration('project_process')

        create_output_links(
            self.job_dir,
            output_cfg,
            dir_with_linked_files
        )

        return output_data_and_archive(
            dir_with_linked_files,
            os.path.join(cfg['project']['input_dir'], self.dataset.name)
        )
