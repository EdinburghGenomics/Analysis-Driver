import os
from egcg_core import executor, util
from analysis_driver import segmentation
from analysis_driver.util import bash_commands
from analysis_driver.pipelines import Pipeline
from analysis_driver.pipelines.common import Cleanup
from analysis_driver.pipelines.qc_gatk4 import GATK4Stage
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import relatedness
from analysis_driver.transfer_data import output_data_and_archive, create_output_links
from analysis_driver.tool_versioning import toolset


class Project(Pipeline):
    """Original GATK3 trio check."""

    toolset_type = 'project_processing'

    def build(self):
        gvcf_files = self.dataset.get_processed_gvcfs()

        if len(gvcf_files) < 2:
            # No need to run as there are not enough gvcfs to process
            return self.stage(Cleanup)
        else:
            return self.trio_check(gvcf_files)

    def trio_check(self, gvcfs):
        sample_ids = [sample['sample_id'] for sample in self.dataset.samples_processed]
        genotype_gvcfs = self.stage(relatedness.GenotypeGVCFs, gVCFs=gvcfs)
        rel = self.stage(relatedness.Relatedness, previous_stages=[genotype_gvcfs])
        peddy = self.stage(relatedness.Peddy, ids=sample_ids, previous_stages=[genotype_gvcfs])
        parse = self.stage(relatedness.ParseRelatedness, ids=sample_ids, parse_method='parse_both', previous_stages=[rel, peddy])
        output = self.stage(Output, previous_stages=[parse])
        cleanup = self.stage(Cleanup, previous_stages=[output])
        return cleanup


class GATK4Project(Project):
    def trio_check(self, gvcfs):
        sample_ids = [sample['sample_id'] for sample in self.dataset.samples_processed]
        db_import = self.stage(GenomicsDBImport, gvcfs=gvcfs)
        genotype_gvcfs = self.stage(GATK4GenotypeGCVFs, previous_stages=[db_import])
        gather = self.stage(GatherVCFs, previous_stages=[genotype_gvcfs])
        rel = self.stage(relatedness.Relatedness, previous_stages=[gather])
        peddy = self.stage(relatedness.Peddy, ids=sample_ids, previous_stages=[gather])
        parse = self.stage(relatedness.ParseRelatedness, ids=sample_ids, parse_method='parse_both',
                           previous_stages=[rel, peddy])
        output = self.stage(Output, previous_stages=[parse])
        cleanup = self.stage(Cleanup, previous_stages=[output])
        return cleanup


class GATK4TrioStage(relatedness.RelatednessStage, GATK4Stage):
    pass


class GenomicsDBImport(GATK4TrioStage):
    gvcfs = segmentation.ListParameter()

    def _run(self):
        cmds = []

        for chunk, region_file in enumerate(self.pipeline.chunk_handler.split_genome_files(os.path.join(self.job_dir, 'genomicsdbimport')).items()):
            genomicsdb = region_file.replace('.bed', '')

            import_cmd = self.gatk_cmd(
                'GenomicsDBImport',
                memory=self.memory,
                reference=None,
                ext='--genomicsdb-workspace-path {db} -L {bed} {inputs}'.format(
                    db=genomicsdb, bed=region_file, inputs=' '.join('-V ' + v for v in self.gvcfs)
                )
            )
            cmds.append(import_cmd)

        return executor.execute(
            *cmds,
            job_name='genomicsdb_import',
            working_dir=self.job_dir,
            cpus=12,
            mem=self.memory
        ).join()


class GATK4GenotypeGCVFs(GATK4TrioStage):
    def _run(self):
        cmds = [
            self.gatk_cmd(
                'GenotypeGVCFs',
                'genotyped_%s.vcf' % c,
                self.memory,
                ext='-V gendb://%s' % os.path.join(self.job_dir, c)
            )
            for c in self.pipeline.chunk_handler.genome_chunks
        ]

        return executor.execute(
            *cmds,
            job_name='genotype_gvcfs',
            working_dir=self.job_dir,
            cpus=12,
            mem=self.memory
        ).join()


class GatherVCFs(GATK4TrioStage):
    def _run(self):
        combine_cmd = self.gatk_picard_cmd(
            'GatherVcfs',
            self.gatk_outfile,
            ext=' '.join('-I genotyped_%s.vcf' % c for c in self.pipeline.chunk_handler.genome_chunks)
        )

        return executor.execute(
            combine_cmd,
            job_name='genotype_cmds',
            working_dir=self.job_dir,
            cpus=12,
            mem=self.memory + 2
        )


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
