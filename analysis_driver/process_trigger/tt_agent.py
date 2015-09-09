__author__ = 'mwham'
import os.path
from time import sleep
from analysis_driver.executor import StreamExecutor
from analysis_driver.config import default as cfg
from analysis_driver.app_logging import get_logger
from .dataset_scanner import is_ready, is_not_ready

app_logger = get_logger('tt_agent')


def rsync(dataset):
    app_logger.info('Starting transfer')
    while is_not_ready(dataset):
        executor = StreamExecutor(
            [
                'rsync',
                '-aqu',
                '--size-only',
                '--partial',
                os.path.join(cfg['raw_dir'], dataset),
                cfg['input_data']
            ]
        )
        executor.run()
        sleep(int(cfg['tt_agent_delay']))

    assert is_ready(dataset)
    app_logger.info('Transfer complete')
