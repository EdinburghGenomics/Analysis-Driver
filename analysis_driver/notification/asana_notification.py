import asana
from .notification_center import Notification


class AsanaNotification(Notification):
    def __init__(self, dataset, config):
        super().__init__(dataset)
        self.client = asana.Client.access_token(config['access_token'])
        self.workspace_id = config['workspace_id']
        self.project_id = config['project_id']
        self._task = None
        self.task_template = {
            'name': self.dataset.name,
            'notes': 'This task has been created because a pipeline run for this dataset had a non-zero exit '
                     'status. Assign this task and close it once the run has been fixed.',
            'projects': [self.project_id]
        }

    @property
    def task(self):
        if self._task is None:
            tasks = list(self.client.tasks.find_all(project=self.project_id, completed=False))
            task_ent = self._get_entity(tasks, self.dataset.name)
            if task_ent is None:
                self.debug('Creating new Asana task')
                task_ent = self._create_task()
            self._task = self.client.tasks.find_by_id(task_ent['id'])
        return self._task

    def end_pipeline(self, exit_status=0):
        if exit_status != 0:
            self._add_comment('Pipeline finished with exit status ' + str(exit_status))

    def crash_report(self, crash_report):
        self._add_comment('Crash report:\n\n' + crash_report)

    def _add_comment(self, text):
        self.client.tasks.add_comment(self.task['id'], text=text)

    def _get_entity(self, collection, name):
        for e in collection:
            if e['name'] == name:
                return e
        self.warning('Task not found: %s' % name)

    def _create_task(self):
        return self.client.tasks.create_in_workspace(self.workspace_id, self.task_template)


if __name__ == '__main__':
    from analysis_driver.dataset_scanner import NoCommunicationDataset
    from analysis_driver.config import default as cfg
    n = AsanaNotification(NoCommunicationDataset('a_test_run'), cfg['notification']['asana'])
    n.end_pipeline(1337)
    n.crash_report('Some\n  stacktrace\n  lines\n')
