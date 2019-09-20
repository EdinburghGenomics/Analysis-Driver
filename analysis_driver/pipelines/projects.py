import os
from egcg_core import executor, util
from analysis_driver import segmentation
from analysis_driver.util import bash_commands
from analysis_driver.pipelines import Pipeline, common
from analysis_driver.pipelines.qc_gatk4 import GATK4Stage, ChunkHandler
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
            return self.stage(common.Cleanup)
        else:
            return self.trio_check(gvcf_files)

    def trio_check(self, gvcfs):
        sample_ids = [sample['sample_id'] for sample in self.dataset.samples_processed]
        genotype_gvcfs = self.stage(relatedness.GenotypeGVCFs, gVCFs=gvcfs)
        rel = self.stage(relatedness.Relatedness, previous_stages=[genotype_gvcfs])
        peddy = self.stage(relatedness.Peddy, ids=sample_ids, previous_stages=[genotype_gvcfs])
        parse = self.stage(relatedness.ParseRelatedness, ids=sample_ids, parse_method='parse_both', previous_stages=[rel, peddy])
        output = self.stage(Output, previous_stages=[parse])
        cleanup = self.stage(common.Cleanup, previous_stages=[output])
        return cleanup


class GATK4Project(Project):
    def __init__(self, dataset):
        super().__init__(dataset)
        self.chunk_handler = ChunkHandler(self.dataset)

    def trio_check(self, gvcfs):
        sample_ids = [sample['sample_id'] for sample in self.dataset.samples_processed]
        db_import = self.stage(GenomicsDBImport, gvcfs=gvcfs)
        genotype_gvcfs = self.stage(GATK4GenotypeGVCFs, stage_name='genotypegvcfs', previous_stages=[db_import])
        gather = self.stage(GatherVCFs, previous_stages=[genotype_gvcfs])
        rel = self.stage(relatedness.Relatedness, previous_stages=[gather])
        peddy = self.stage(relatedness.Peddy, ids=sample_ids, previous_stages=[gather])
        parse = self.stage(relatedness.ParseRelatedness, ids=sample_ids, parse_method='parse_both',
                           previous_stages=[rel, peddy])
        output = self.stage(Output, previous_stages=[parse])
        cleanup = self.stage(common.Cleanup, previous_stages=[output])
        return cleanup


class GATK4TrioStage(relatedness.RelatednessStage, GATK4Stage):
    def chunk_file(self, chunk, suffix=''):
        return os.path.join(self.job_dir, self.dataset.name + '_region_%s-%s-%s' % chunk + suffix)


class GenomicsDBImport(GATK4TrioStage):
    gvcfs = segmentation.ListParameter()

    def _run(self):
        cmds = []

        for chunk, region_file in self.pipeline.chunk_handler.split_genome_files(self.job_dir).items():
            import_cmd = self.gatk_cmd(
                'GenomicsDBImport',
                memory=self.memory(self.gvcfs),
                reference=None,
                ext='--genomicsdb-workspace-path {db} -L {bed} {inputs}'.format(
                    db=self.chunk_file(chunk), bed=region_file, inputs=' '.join('-V ' + v for v in self.gvcfs)
                )
            )
            cmds.append(import_cmd)

        return executor.execute(
            *cmds,
            job_name='genomicsdb_import',
            working_dir=self.job_dir,
            cpus=12,
            mem=self.memory(self.gvcfs)
        ).join()


class GATK4GenotypeGVCFs(GATK4TrioStage):
    def _run(self):
        cmds = [
            self.gatk_cmd(
                'GenotypeGVCFs',
                'genotyped_%s-%s-%s.vcf' % c[0],
                12,
                ext='-V gendb://%s' % self.chunk_file(c[0])
            )
            for c in self.pipeline.chunk_handler.genome_chunks
        ]

        return executor.execute(
            *cmds,
            job_name='genotype_gvcfs',
            working_dir=self.job_dir,
            cpus=2,
            mem=12
        ).join()


class GatherVCFs(GATK4TrioStage):
    def _run(self):
        combined_vcf = os.path.join(self.job_dir, 'combined.vcf')
        combine_cmd = self.gatk_picard_cmd(
            'GatherVcfs',
            combined_vcf,
            reference=None,
            ext=' '.join('-I genotyped_%s-%s-%s.vcf' % c[0] for c in self.pipeline.chunk_handler.genome_chunks)
        )

        gather = executor.execute(
            combine_cmd,
            job_name='gather_vcfs',
            working_dir=self.job_dir,
            cpus=12,
            mem=self.memory(self.pipeline.chunk_handler.genome_chunks)
        ).join()

        reorder_cmd = self.gatk_picard_cmd(
            'SortVcf',
            self.gatk_outfile,
            reference=None,
            input=combined_vcf,
            ext='-SD ' + os.path.splitext(self.dataset.reference_genome)[0] + '.dict'
        )

        reorder = executor.execute(
            reorder_cmd,
            job_name='reorder_vcf',
            working_dir=self.job_dir,
            cpus=12,
            mem=self.memory(self.pipeline.chunk_handler.genome_chunks)
        ).join()

        return gather + reorder


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
