from analysis_driver import rest_communication


class SampleReview():

    def __init__(self, name):
        self.name = name
        self.run_elements = self._read_data()

    def _read_data(self):
        return rest_communication.get_documents(
                'run_elements',
                where={'sample_id': self.name}
            )

    def print_data(self):
        for r in self.run_elements:
            print(r)