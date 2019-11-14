import pytest
from os.path import join
from unittest.mock import Mock, patch
from analysis_driver import segmentation
from analysis_driver.dataset import RunDataset
from analysis_driver.exceptions import SequencingRunError, PipelineError
from analysis_driver.pipelines import pipeline, bcbio, demultiplexing, human_variant_calling_gatk4, projects, qc,\
    qc_gatk4, rapid, variant_calling, variant_calling_gatk4
from tests.test_analysisdriver import TestAnalysisDriver, NamedMock


class ExceptionStage(segmentation.Stage):
    def _run(self):
        raise SequencingRunError('Aborted')


class SuccessfulStage(segmentation.Stage):
    def _run(self):
        pass


class FailingStage(segmentation.Stage):
    def _run(self):
        return 1


class TestPipeline(TestAnalysisDriver):
    @staticmethod
    def _pipeline_with_stage(stage_class):
        patches = []
        for function_to_patch in ['start', 'get_stage', 'start_stage', 'end_stage', 'resolve_pipeline_and_toolset']:
            p = patch.object(RunDataset, function_to_patch)
            p.start()
            patches.append(p)

        dataset = RunDataset(name='a_sample')
        patched_pipeline = patch.object(
            RunDataset,
            'pipeline',
            new=Mock(build=Mock(return_value=stage_class(dataset=dataset, pipeline=Mock())))
        )
        patched_pipeline.start()
        patches.append(patched_pipeline)

        try:
            return pipeline(dataset)
        finally:
            for p in patches:
                p.stop()

    def test_pipeline_with_exceptions(self):
        with pytest.raises(SequencingRunError):
            self._pipeline_with_stage(ExceptionStage)

    def test_successful_pipeline(self):
        assert self._pipeline_with_stage(SuccessfulStage) == 0

    def test_failing_pipeline(self):
        with pytest.raises(PipelineError):
            self._pipeline_with_stage(FailingStage)

    def test_pipelines(self):
        sample_pipelines = (
            bcbio.BCBioVarCalling, human_variant_calling_gatk4.HumanVarCallingGATK4, qc.QC, qc_gatk4.QCGATK4,
            variant_calling.VarCalling, variant_calling_gatk4.VarCallingGATK4
        )

        run_dataset = NamedMock(real_name='test_run', type='run')
        demultiplexing.Demultiplexing(run_dataset).build()
        rapid.Rapid(run_dataset).build(Mock())

        sample_dataset = NamedMock(
            real_name='test_sample',
            user_sample_id='uid',
            type='sample',
            reference_genome=join(self.assets_path, 'genome.fa'),
            genome_dict={'snpEff': 'path/to/snpEff'}
        )
        for cls in sample_pipelines:
            cls(sample_dataset).build()

        project_dataset = NamedMock(
            real_name='test_project',
            type='project',
            reference_genome=join(self.assets_path, 'genome.fa'),
            get_processed_gvcfs=Mock(return_value=['gvcf1', 'gvcf2']),
            samples_processed=[{'sample_id': 'sample_1'}, {'sample_id': 'sample_2'}]
        )
        projects.Project(project_dataset).build()
