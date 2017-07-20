from egcg_core import notifications, clarity
from pyclarity_lims.lims import Step
from analysis_driver.exceptions import PipelineError


class LimsNotification():
    def __init__(self, name):
        self.name = name
        self.sample_list = clarity.get_list_of_samples(self.name)
        self.artifacts = [sample.artifact for sample in self.sample_list]
        self.step = None

    @property
    def lims(self):
        return clarity.connection()

    @property
    def postseqlab_stages(self):
        return self.lims.get_workflows(name='PostSeqLab EG 1.0 WF')[0].stages

    def get_stage(self, stage_name):
        stage = [stage for stage in self.postseqlab_stages if stage == stage_name]
        if not len(stage) == 1:
            raise PipelineError('Can not find %s from Workflow Stages' % stage_name)
        return stage[0]

    def route_sample_to_data_processing(self):
        self.lims.route_artifacts(self.artifacts, stage_uri=self.get_stage('Data Processing EG 1.0 ST').uri)

    def create_step(self):
        self.step = Step.create(self.lims, protocol_step=self.get_stage('Data Processing EG 1.0 ST').step, inputs=[self.artifacts], container_type_name='Tube')
        self.step.advance()

    def assign_next_and_advance_step(self):
        assert self.step.current_state == 'Assign Next Steps'
        next_step_uri = self.get_stage('Sample Review EG 1.0 ST').uri
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

    def end_sample_pipeline(self):
        self.assign_next_and_advance_step()

    def fail_and_reset_sample_pipeline(self):
        self.remove_sample_from_workflow()


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
