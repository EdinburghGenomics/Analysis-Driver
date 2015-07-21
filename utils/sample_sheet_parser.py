import os.path
import csv


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
        last_sample = None
        for name, sample_project in self.sample_projects.items():
            last_sample = None
            for sample in sample_project.samples:
                try:
                    if len(sample.barcode) == len(last_sample.barcode):
                        pass
                    else:
                        raise ValueError(
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
            # 'workspace',
            # 'EdGen_Analysis_Driver',
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
