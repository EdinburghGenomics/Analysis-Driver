import os
from egcg_core import executor
from egcg_core.util import find_file, query_dict
from analysis_driver import segmentation
from analysis_driver.util import bash_commands
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.config import default as cfg, output_file_config
from analysis_driver.quality_control import Relatedness, Peddy, GenotypeGVCFs, ParseRelatedness
from analysis_driver.exceptions import PipelineError
from analysis_driver.transfer_data import output_data_and_archive, create_output_links
from analysis_driver.tool_versioning import toolset

toolset_type = 'project_processing'
name = 'project'


def build_pipeline(dataset):
    sample_ids = [sample['sample_id'] for sample in dataset.samples_processed]
    project_source = os.path.join(cfg.query('project', 'input_dir'), dataset.name)
    gvcf_files = []
    for sample in dataset.samples_processed:
        # Only check if we have gvcf when the samples have been through human processing that generate a gvcf
        if query_dict(sample, 'aggregated.most_recent_proc.pipeline_used.name') == 'bcbio':
            gvcf_file = find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
            if not gvcf_file:
                raise PipelineError('Unable to find gVCF file for sample %s in %s' % (sample['sample_id'], project_source))
            gvcf_files.append(gvcf_file)

    if len(gvcf_files) < 2:
        # No need to run as there are not enough gvcfs to process
        cleanup = Cleanup(dataset=dataset)
    else:
        genotype_gvcfs = GenotypeGVCFs(dataset=dataset, gVCFs=gvcf_files)
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
        os.makedirs(dir_with_linked_files, exist_ok=True)
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
        output_file_config.set_pipeline_type('project_process')
        toolset.write_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))

        create_output_links(
            self.job_dir,
            output_file_config,
            dir_with_linked_files,
            project_id=self.dataset.name
        )

        return output_data_and_archive(
            dir_with_linked_files,
            os.path.join(cfg.query('project', 'output_dir'), self.dataset.name)
        )
