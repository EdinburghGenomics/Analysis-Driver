from egcg_core.exceptions import EGCGError


class AnalysisDriverError(EGCGError):
    pass


class PipelineError(AnalysisDriverError):
    pass


class SequencingRunError(AnalysisDriverError):
    pass

