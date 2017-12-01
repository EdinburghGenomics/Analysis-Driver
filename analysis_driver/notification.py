from egcg_core import notifications, clarity
from pyclarity_lims.entities import Step


class LimsNotification:
    def __init__(self, name):
        self.name = name
        self.sample = clarity.get_sample(self.name)
        self.artifact = self.sample.artifact
        self.data_processing_stage = clarity.get_workflow_stage('PostSeqLab EG 1.0 WF', stage_name='Data Processing EG 1.0 ST')
        self.sample_review_stage = clarity.get_workflow_stage('PostSeqLab EG 1.0 WF', stage_name='Sample Review EG 1.0 ST')
        self.step = None

    @property
    def lims(self):
        return clarity.connection()

    def route_sample_to_data_processing(self):
        self.lims.route_artifacts([self.artifact], stage_uri=self.data_processing_stage.uri)

    def create_step(self):
        self.step = Step.create(self.lims, protocol_step=self.data_processing_stage.step, inputs=[self.artifact])
        self.step.advance()

    def assign_next_and_advance_step(self):
        assert self.step.current_state == 'Assign Next Steps'
        next_step_uri = self.sample_review_stage.uri
        self.step.actions.next_actions[0]['step-uri'] = next_step_uri
        self.step.actions.next_actions[0]['action'] = 'nextstep'
        self.step.actions.put()
        self.step.advance()

    def remove_sample_from_workflow(self):
        assert self.step.current_state == 'Assign Next Steps'
        self.step.actions.next_actions[0]['action'] = 'remove'
        self.step.actions.put()
        self.step.advance()

    def start_sample_pipeline(self):
        self.route_sample_to_data_processing()
        self.create_step()


class NotificationCentre(notifications.NotificationCentre):
    def start_pipeline(self):
        self.notify('Started pipeline', ('log', 'email'))

    def start_stage(self, stage_name):
        self.notify('Started stage ' + stage_name, ('log',))

    def end_stage(self, stage_name, exit_status=0):
        subs = ['log']
        if exit_status:
            subs.append('email')
        self.notify("Stage '%s' finished with exit status %s" % (stage_name, exit_status), subs)

    def end_pipeline(self, exit_status):
        subs = ['log', 'email']
        if exit_status:
            subs.append('asana')
        self.notify('Pipeline finished with exit status ' + str(exit_status), subs)

    def crash_report(self, stacktrace):
        self.notify_all(stacktrace)
