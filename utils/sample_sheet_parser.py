import os.path
import csv


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
        d.setdefault(sample_id, []).append([lane, sample_name, first_time_pos, project_name])

    return d


def get_sample_project(d):
    #  the first one for instance
    return d.values()[0][0][3]


if __name__ == '__main__':
    print(read_sample_sheet('/home/U008/mwham/input_data/150424_E00307_0017_AH3KGTCCXX'))

