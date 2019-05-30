import os
from egcg_core import executor, util
from analysis_driver import segmentation
from analysis_driver.util import bash_commands
from analysis_driver.pipelines import Pipeline
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import Relatedness, Peddy, GenotypeGVCFs, ParseRelatedness
from analysis_driver.exceptions import PipelineError
from analysis_driver.transfer_data import output_data_and_archive, create_output_links
from analysis_driver.tool_versioning import toolset


class Project(Pipeline):
    toolset_type = 'project_processing'

    def build(self):
        sample_ids = [sample['sample_id'] for sample in self.dataset.samples_processed]
        project_source = os.path.join(cfg.query('project', 'input_dir'), self.dataset.name)
        gvcf_files = []
        for sample in self.dataset.samples_processed:
            # Only check if we have gvcf when the samples have been through human processing that generate a gvcf
            if util.query_dict(sample, 'aggregated.most_recent_proc.pipeline_used.name') == 'bcbio':
                gvcf_file = util.find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
                if not gvcf_file:
                    raise PipelineError('Unable to find gVCF file for sample %s in %s' % (sample['sample_id'], project_source))
                gvcf_files.append(gvcf_file)

        if len(gvcf_files) < 2:
            # No need to run as there are not enough gvcfs to process
            cleanup = self.stage(Cleanup)
        else:
            genotype_gvcfs = self.stage(GenotypeGVCFs, gVCFs=gvcf_files)
            relatedness = self.stage(Relatedness, previous_stages=[genotype_gvcfs])
            peddy = self.stage(Peddy, ids=sample_ids, previous_stages=[genotype_gvcfs])
            parse = self.stage(ParseRelatedness, ids=sample_ids, parse_method='parse_both', previous_stages=[relatedness, peddy])
            output = self.stage(Output, previous_stages=[parse])
            cleanup = self.stage(Cleanup, previous_stages=[output])
        return cleanup


class Output(segmentation.Stage):
    def _run(self):
        dir_with_linked_files = os.path.join(self.job_dir, 'relatedness_outfiles')
        os.makedirs(dir_with_linked_files, exist_ok=True)
        toolset.write_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))

        create_output_links(
            self.job_dir,
            'project_process',
            dir_with_linked_files,
            project_id=self.dataset.name
        )

        md5_exit_status = executor.execute(
            *[bash_commands.md5sum(f) for f in sorted(util.find_files(dir_with_linked_files, '*'))],
            job_name='md5sum',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()

        return md5_exit_status + output_data_and_archive(
            dir_with_linked_files,
            os.path.join(cfg.query('project', 'output_dir'), self.dataset.name)
        )
