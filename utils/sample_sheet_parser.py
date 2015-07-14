import os.path
import csv
#from util.logger import AppLogger


def read_sample_sheet(file_path):
    """
    Read sample IDs and sample names from SampleSheet.csv
    :param file_path: Full path to SampleSheet.csv
    :return: A dict of sample IDs and sample names
    """

    d = {}
    tmp_list = []
    filename = os.path.join(file_path, 'SampleSheet.csv')  # hardcoded name
    reader = csv.reader(open(filename, 'r'))

    # read lines until [Data] marker
    while not next(reader)[0].startswith('[Data]'):
        pass
    next(reader)  # one more to ignore header

    # store the rest of the csv file in a list
    for row in reader:
        # store only non empty lines
        if any(row):
            tmp_list.append(row)
    # get the position in the CSV file (required for file naming)
    first_time = {}
    for row in tmp_list:
        if row[1] not in first_time:
            first_time[row[1]] = row[0]

    # store keys (sampleID) and values (Lane,SampleName, firstTimePos, projectName)
    # {
    #     '10015AT0001L01': [
    #         ['1', '50293', '1', '10015AT']
    #     ],
    #     '10015AT0002L01': [
    #         ['2', '128864', '2', '10015AT']
    #     ],
    #     '10015AT0003L01': [
    #         ['3', '172937', '3', '10015AT']
    #     ],
    #     '10015AT0004L01': [
    #         ['4', 'Na12878', '4', '10015AT'],
    #         ['6', 'Na12878', '4', '10015AT'],
    #         ['7', 'Na12878', '4', '10015AT']
    #     ],
    #     '10015AT0005L01': [
    #         ['5', 'PhiX', '5', '10015AT'],
    #         ['8', 'PhiX', '5', '10015AT']
    #     ]
    # }

    for row in tmp_list:
        sample_id = row[1]
        lane = row[0]
        sample_name = row[2]
        first_time_pos = first_time[row[1]]
        project_name = row[7]
        # project_name = row[5]
        d.setdefault(sample_id, []).append([lane, sample_name, first_time_pos, project_name])

    return d


def get_sample_project(d):
    #  the first one for instance
    return list(d.values())[0][0][3]


class SampleSheet:
    def __init__(self, data_dir):
        self.file = open(os.path.join(data_dir, 'SampleSheet.csv'), 'r')
        # read lines until [Data] marker
        while not next(self.file).startswith('[Data]'):
            pass

        self.sample_projects = {}  # {name: samples} {str: Sample}
        self._populate(self.sample_projects)

    def _populate(self, samples):
        reader = csv.DictReader(self.file)
        for row in reader:
            if any(row):
                sample_project = row['Sample_Project']
                new_sample = Sample(
                    sample_project,
                    row['Lane'],
                    row['Sample_ID'],
                    row['Sample_Name'],
                    row['I7_Index_ID'],
                    row['index'],
                    row['Sample_Plate'],
                    row['Sample_Well'],
                    row['Description']
                )
                try:
                    samples[sample_project].add_sample(new_sample)
                except KeyError:
                    samples[sample_project] = SampleProject(sample_project, new_sample)

    def check_barcodes(self):
        errors = []
        for name, sample_project in self.sample_projects.items():
            last_sample = None
            for sample in sample_project.samples:
                try:
                    if len(sample.barcode) == len(last_sample.barcode):
                        pass
                    else:
                        raise  ValueError(
                            'Odd barcode length for %s: %s (%s) in sample project \'%s\' ' % (
                                sample.id, sample.barcode, len(sample.barcode), name
                            )
                        )
                except AttributeError:
                    pass
                finally:
                    last_sample = sample

        return len(last_sample.barcode)


class SampleProject:
    def __init__(self, name, new_sample):
        self.name = name
        self.samples = [new_sample]

    def add_sample(self, sample):
        if sample.sample_project != self.name:
            raise AssertionError
        else:
            self.samples.append(sample)


class Sample:
    def __init__(self, sample_project, lane, sample_id, sample_name, i7_index_id,
                 index=None, sample_plate=None, sample_well=None, description=None):
        self.sample_project = sample_project
        self.lane = lane
        self.id = sample_id
        self.name = sample_name
        self.index_id = i7_index_id
        self.barcode = index
        self.plate = sample_plate
        self.well = sample_well
        self.description = description



if __name__ == '__main__':
    import os
    sheet = SampleSheet(
        os.path.join(
            os.getenv('HOME'),
#            'workspace',
#            'EdGen_Analysis_Driver',
            'input_data',
            '150424_E00307_0017_AH3KGTCCXX'
        )
    )
    print(
        [
            (
                name, [sample.id for sample in sample_project.samples]
            ) for name, sample_project in sheet.sample_projects.items()
        ]
    )
    print(sheet.check_barcodes())
