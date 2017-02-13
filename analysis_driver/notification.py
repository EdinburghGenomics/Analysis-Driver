from egcg_core import notifications


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
