import unittest
from unittest.mock import Mock, patch

import luigi
import pytest
from luigi import Task, Event

from analysis_driver import segmentation
from analysis_driver.dataset import RunDataset
from analysis_driver.exceptions import SequencingRunError, PipelineError
from analysis_driver.pipelines import pipeline
from tests import TestAnalysisDriver
from tests.test_analysisdriver import NamedMock


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

    def _pipeline_with_stage(self, stage_class):
        dataset = RunDataset(name='a_sample')
        patches = []
        for function_to_patch in ['start', 'get_stage', 'start_stage', 'end_stage', 'resolve_pipeline_and_toolset']:
            p = patch.object(RunDataset, function_to_patch)
            p.start()
            patches.append(p)

        dataset.pipeline = Mock(build_pipeline=Mock(return_value=stage_class(dataset=dataset)))
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

    def test_failling_pipeline(self):
        with pytest.raises(PipelineError):
            self._pipeline_with_stage(FailingStage)


