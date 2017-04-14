import os
import yaml
from egcg_core import executor, clarity, util
from analysis_driver import segmentation
from analysis_driver import quality_control as qc
from analysis_driver.exceptions import PipelineError
from analysis_driver.pipelines.common import bcbio_prepare_sample, link_results_files, output_data, cleanup
from analysis_driver.util import bash_commands
from analysis_driver.config import default as cfg
from analysis_driver.reader.version_reader import write_versions_to_yaml
from analysis_driver.transfer_data import prepare_sample_data


class BCBioStage(segmentation.Stage):
    @property
    def fastq_pair(self):
        return util.find_files(self.job_dir, 'merged', self.dataset.user_sample_id + '_R?.fastq.gz')


class MergeFastqs(BCBioStage):
    def _run(self):
        fastq_files = prepare_sample_data(self.dataset)
        bcbio_prepare_sample(self.job_dir, self.dataset.name, fastq_files)
        return 0


class FastQC(BCBioStage):
    previous_stages = [MergeFastqs]

    def _run(self):
        return executor.execute(
            *[bash_commands.fastqc(fastq_file) for fastq_file in self.fastq_pair],
            job_name='fastqc2',
            working_dir=self.job_dir,
            cpus=1,
            mem=2
        ).join()


class BCBioAndQC(segmentation.BasicStage):
    previous_stages = [
        (FastQC, {'previous_stages': [FastQC]}),
        (qc.ContaminationCheck, {'previous_stages': [FastQC]}),
        (qc.ContaminationBlast, {'previous_stages': [FastQC]}),
        (qc.GenotypeValidation, {'previous_stages': [FastQC]}),
        (BCBio, {'previous_stages': [FastQC]})
    ]


class PostBCBioQC(segmentation.BasicStage):
    @property
    def previous_stages(self):
        return [
            (qc.GenderValidation, {'vcf_file': self.vcf_file, 'previous_stages': BCBioAndQC}),
            (qc.VCFStats, {'vcf_file': self.vcf_file, 'previous_stages': BCBioAndQC}),
            (qc.VerifyBamID, {'bam_file': self.bam_file, 'previous_stages': BCBioAndQC}),
            (qc.SamtoolsDepth, {'bam_file': self.bam_file, 'previous_stages': BCBioAndQC})
        ]

    @property
    def vcf(self):
        return util.find_file(
            self.job_dir,
            'samples_%s-merged' % self.dataset.name,
            'final',
            '*_%s' % self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-joint-gatk-haplotype-joint.vcf.gz'
        )

    @property
    def bam(self):
        return util.find_file(
            'samples_%s-merged' % self.dataset.name,
            'final',
            self.dataset.user_sample_id,
            self.dataset.user_sample_id + '-ready.bam'
        )


class Output(BCBioStage):
    previous_stages = [PostBCBioQC]

    def _run(self):
        # link the bcbio file into the final directory
        dir_with_linked_files = link_results_files(self.dataset.name, self.job_dir, 'bcbio')
        write_versions_to_yaml(os.path.join(dir_with_linked_files, 'program_versions.yaml'))
        return output_data(self.dataset, self.job_dir, self.dataset.name, dir_with_linked_files)


class Cleanup(BCBioStage):
    previous_stages = [Output]

    def _run(self):
        return cleanup(self.dataset.name)


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
            '..', '..', 'etc', 'bcbio_alignment_%s_%s.yaml' % (self.genome_version, analysis_type)
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
            os.path.join(cfg['tools']['bcbio'], 'bin', 'bcbio_nextgen.py'),
            '-w template',
            run_template,
            bcbio_dir,
            bcbio_dir + '.csv'
        ] + self.fastq_pair

        run_yaml = os.path.join(bcbio_dir, 'config', 'samples_' + self.dataset.name + '-merged.yaml')
        bcbio_cmd = bash_commands.bcbio(run_yaml, os.path.join(bcbio_dir, 'work'), threads=16)

        prep_status = executor.execute(' '.join(sample_prep), env='local').join()
        self.info('BCBio sample prep exit status: ' + str(prep_status))

        ext_sample_id = clarity.get_user_sample_name(self.dataset.name, lenient=True)

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
