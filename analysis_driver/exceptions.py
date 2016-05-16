class AnalysisDriverError(Exception):
    pass


class PipelineError(AnalysisDriverError):
    pass


class SequencingRunError(AnalysisDriverError):
    pass


class RestCommunicationError(AnalysisDriverError):
    pass


class LimsCommunicationError(AnalysisDriverError):
    pass
