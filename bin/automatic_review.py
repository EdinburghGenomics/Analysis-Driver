from analysis_driver import rest_communication


class SampleReview():

    def __init__(self, name):
        self.name = name
        self.run_elements = self.read_data()

    def read_data(self):
        return rest_communication.get_documents(
            'aggregate/samples',
            depaginate=True,
            max_results=100,
            embedded={"run_elements":1},
            match={"sample_id":self.name}
        )

    def print_data(self):
        for r in self.run_elements:
            print(r)

s = SampleReview("10015AT0009")
s.print_data()