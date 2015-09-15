__author__ = 'mwham'
import os
from analysis_driver.config import default as cfg
from analysis_driver.dataset_scanner import dataset_status, switch_status
from . import tt_agent


def trigger(dataset, use_intermediate_dir):
    if use_intermediate_dir:
        assert dataset_status(dataset) in ('new', 'new, rta complete')
        switch_status(dataset, 'transferring')
        tt_agent.transfer_to_int_dir(
            dataset,
            cfg['input_dir'],
            cfg['intermediate_dir'],
            repeat_delay=int(cfg['tt_agent_delay'])
        )
        dataset_dir = cfg['intermediate_dir']
    else:
        assert dataset_status(dataset) == 'new, rta_complete'
        dataset_dir = cfg['input_dir']

    switch_status(dataset, 'active')

    from analysis_driver import driver
    driver.pipeline(os.path.join(dataset_dir, dataset))

    switch_status(dataset, 'complete')
