import os
import yaml
from egcg_core import executor, clarity, util
from analysis_driver import segmentation, quality_control as qc
from analysis_driver.exceptions import PipelineError
from analysis_driver.pipelines import common
from analysis_driver.util import bash_commands
from analysis_driver.tool_versioning import toolset

toolset_type = 'human_sample_processing'
name = 'bcbio'


class BCBioStage(segmentation.Stage):
    @property
    def fastq_pair(self):
        return os.path.join(self.job_dir, 'merged', self.dataset.user_sample_id + '_R?.fastq.gz')

    @property
    def vcf_path(self):
        return os.path.join(
            self.job_dir,
            'samples_%s-merged' % self.dataset.name,
            'final',
            '*_%s' % self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-joint-gatk-haplotype-joint.vcf.gz'
        )

    @property
    def bam_path(self):
        return os.path.join(
            self.job_dir,
            'samples_%s-merged' % self.dataset.name,
            'final',
            self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-ready.bam'
        )


class BCBio(BCBioStage):
    def _run(self):
        analysis_type = clarity.get_sample(self.dataset.name).udf.get('Analysis Type')

        if not analysis_type:
            analysis_type = 'gatk'
        elif analysis_type.endswith('gatk'):
            analysis_type = 'gatk'
        elif analysis_type.endswith('freebayes'):
            analysis_type = 'freebayes'
        else:
            self.error('Unknown Analysis type %s' % analysis_type)
            return 1

        run_template = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            '..', '..', 'etc', 'bcbio_alignment_%s_%s.yaml' % (self.dataset.genome_version, analysis_type)
        )
        if not os.path.isfile(run_template):
            raise PipelineError(
                'Could not find BCBio run template ' + run_template + '. Is the correct genome set?'
            )

        original_dir = os.getcwd()
        os.chdir(self.job_dir)
        self.debug('Using merged fastqs: %s', self.fastq_pair)

        bcbio_dir = os.path.join(self.job_dir, 'samples_' + self.dataset.name + '-merged')

        sample_prep = [
            os.path.join(toolset['bcbio'], 'bin', 'bcbio_nextgen.py'),
            '-w template',
            run_template,
            bcbio_dir,
            bcbio_dir + '.csv'
        ] + util.find_files(self.fastq_pair)

        run_yaml = os.path.join(bcbio_dir, 'config', 'samples_' + self.dataset.name + '-merged.yaml')
        bcbio_cmd = bash_commands.bcbio(run_yaml, os.path.join(bcbio_dir, 'work'), threads=16)

        prep_status = executor.execute(' '.join(sample_prep), env='local').join()
        self.info('BCBio sample prep exit status: ' + str(prep_status))

        ext_sample_id = self.dataset.user_sample_id or self.dataset.name

        with open(run_yaml, 'r') as i:
            run_config = yaml.load(i)
        run_config['fc_name'] = ext_sample_id
        with open(run_yaml, 'w') as o:
            o.write(yaml.safe_dump(run_config, default_flow_style=False))

        bcbio_executor = executor.execute(
            bcbio_cmd,
            prelim_cmds=bash_commands.export_env_vars(),
            job_name='bcb%s' % self.dataset.name,
            working_dir=self.job_dir,
            cpus=12,
            mem=64
        )
        os.chdir(original_dir)
        return bcbio_executor.join()


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    merge_fastqs = stage(common.MergeFastqs)
    fastqc = stage(common.FastQC, previous_stages=[merge_fastqs])
    bcbio = stage(BCBio, previous_stages=[fastqc])

    contam_check = stage(qc.FastqScreen, fq_pattern=bcbio.fastq_pair, previous_stages=[fastqc])
    blast = stage(qc.Blast, fastq_file=bcbio.fastq_pair.replace('?', '1'), previous_stages=[fastqc])
    geno_val = stage(qc.GenotypeValidation, fq_pattern=bcbio.fastq_pair, previous_stages=[fastqc])
    bcbio_and_qc = [bcbio, fastqc, contam_check, blast, geno_val]

    gender_val = stage(qc.GenderValidation, vcf_file=bcbio.vcf_path, previous_stages=bcbio_and_qc),
    vcfstats = stage(qc.VCFStats, vcf_file=bcbio.vcf_path, previous_stages=bcbio_and_qc),
    verify_bam_id = stage(qc.VerifyBamID, bam_file=bcbio.bam_path, previous_stages=bcbio_and_qc),
    samtools_depth = stage(qc.SamtoolsDepth, bam_file=bcbio.bam_path, previous_stages=bcbio_and_qc)
    post_bcbio_qc = [gender_val, vcfstats, verify_bam_id, samtools_depth]

    output = stage(common.SampleDataOutput, previous_stages=post_bcbio_qc, output_fileset='bcbio')
    cleanup = stage(common.Cleanup, previous_stages=[output])

    return cleanup
