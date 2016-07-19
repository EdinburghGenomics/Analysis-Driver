from threading import Thread
from egcg_core.app_logging import AppLogger


class QualityControl(AppLogger, Thread):
    def __init__(self, dataset, working_dir):
        self.exception = None
        self.exit_status = 0
        self.dataset = dataset
        self.working_dir = working_dir
        Thread.__init__(self)
