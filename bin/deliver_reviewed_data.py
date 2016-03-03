import datetime

from analysis_driver import executor, clarity
from analysis_driver.exceptions import AnalysisDriverError

__author__ = 'tcezard'
import sys
import os
import argparse
import logging

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.app_logging import AppLogger


class DataDelivery(AppLogger):

    def __init__(self):
        self.all_commands = []


    def get_deliverable_projects_samples(self):
        project_to_samples = {}
        #Get sample that have been review from REST API but not marked as status release
        return project_to_samples

    def append_create_batch_delivery_folder(self, deliver_dir, project):
        batch = datetime.date.today().isoformat()
        batch_delivery_folder = os.path.join(deliver_dir, project, batch)
        command = 'mkdir -p %s'%batch_delivery_folder
        self.all_commands.append(command)
        return batch_delivery_folder

    def append_create_sample_folder(self, batch_delivery_folder, sample):
        sample_folder = os.path.join(batch_delivery_folder, sample)
        command = 'mkdir -p %s'%sample_folder
        self.all_commands.append(command)
        return sample_folder

    def append_move_file_to_sample_folder(self, file_to_move, sample_folder):
        command = 'mv %s %s'%(file_to_move, sample_folder)
        self.all_commands.append(command)

    def get_list_of_file_to_move(self, sample_name):
        species = clarity.get_species_from_sample(sample_name)


    def deliver_data(self):
        output_dir = cfg.query('output_dir')

        project_to_samples = self.get_deliverable_projects_samples()

        for project in project_to_samples:
            samples = project_to_samples.get(project)
            batch_delivery_folder = self.create_batch_delivery_folder(output_dir, project)
            for sample in samples:
                origin_sample_dir = os.path.join(output_dir, project, sample)
                if not os.path.isdir(origin_sample_dir):
                    raise AnalysisDriverError('Directory for sample %s in project %s does not exist'%(sample, project))
                sample_folder = self.append_create_sample_folder(batch_delivery_folder, sample)
                for file_to_move in list_of_file_to_move:
                    if not os.path.isfile(os.path.join(origin_sample_dir, file_to_move)):
                        raise AnalysisDriverError('File %s for sample %s in project %s does not exist'%(file_to_move, sample, project))
                    self.append_move_file_to_sample_folder(file_to_move, sample_folder)


        if self.dry_run:
            print('\n'.join(self.all_commands))
        else:
            self.mark_samples_as_released()
            for command in self.all_commands:
                exit_status = executor.execute([command], env='local').join()
                if exit_status != 0:
                    raise AnalysisDriverError('command %s exited with status %s'%(command, exit_status))





def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dry_run', action='store_true')
    p.add_argument('--debug', action='store_true')
    p.add_argument('--work_dir', type=str, required=True)
    args = p.parse_args()

    if args.debug:
        log_cfg.default_level = logging.DEBUG
        log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)

    cfg.merge(cfg['sample'])


if __name__ == '__main__':
    main()
