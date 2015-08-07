__author__ = 'mwham'

class ProcessTriggerError(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
