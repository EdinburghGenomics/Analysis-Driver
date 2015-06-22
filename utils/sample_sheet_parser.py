import csv


def read_sample_sheet(file_path):
    """
    Reads the sample ids and sample names
    :param file_path: Full path to SampleSheet.csv
    :return: sample ids and sample names
    :rtype: tuple[int, list[str]]
    """

    d = {}
    csv_data_lines = []
    # open the file with csv (hardcoded name)
    reader = csv.reader(open(file_path + '/SampleSheet.csv', 'rb'))

    # read lines until [Data] marker
    while not next(reader)[0].startswith('[Data]'):
        pass
    next(reader)  # one more - ignore names

    # store the rest of the csv file in a list
    for row in reader:
        if any(row):  # store only non empty lines
            csv_data_lines.append(row)

    # get the position in the CSV file (required for file naming)
    first_time = {}
    for row in csv_data_lines:
        if row[1] not in first_time:
            first_time[row[1]] = row[0]

        # Store keys (sampleID) and values (Lane,SampleName, firstTimePos, projectName)

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

    for row in csv_data_lines:
        d.setdefault(row[1], []).append([row[0], row[2], first_time[row[1]], row[7]])

    return len(csv_data_lines), d


def get_sample_project(d):
    """
    Gives the sample_project identifier for a given sample sheet
    :param d:
    :return:
    """
    # The first one for instance
    return d.values()[0][0][3]
