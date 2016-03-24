class AnalysisDriverError(Exception):
    pass


class PipelineError(AnalysisDriverError):
    pass


class SequencingRunError(AnalysisDriverError):
    pass
