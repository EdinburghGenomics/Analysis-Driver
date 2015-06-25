import csv
# reads the sampleIds and sampleNames 
# Input: Full path to SampleSheet.csv
# Outout: It returns to lists
#         1 - list of sampleIds
#         2-  list of sampleNames

def read_sample_sheet(file_path):
    """
    Read sample IDs and sample names
    :param file_path: Full path to SampleSheet.csv
    :rtype: tuple[int, dict[str[list[list[str, str, str]]]]]
    """

    d = {}
    tmp_list = []
    filename = file_path + '/SampleSheet.csv'  # hardcoded name
    reader = csv.reader(open(filename, 'rb'))

    # read lines until [Data] marker
    while not next(reader)[0].startswith('[Data]'):
        pass
    next(reader)  # ignore names

    # store the rest of the csv file in a list
    for row in reader:
        # store only the non empty lines
        if any(row):
            tmp_list.append(row)
    # get the position in the CSV file (required for file naming)
    first_time = {}
    for row in tmp_list:
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
    #         ['7', 'Na12878', '4', '10015AT']],
    #     '10015AT0005L01': [
    #         ['5', 'PhiX', '5', '10015AT'],
    #         ['8', 'PhiX', '5', '10015AT']
    #     ]
    # }

    for row in tmp_list:
        d.setdefault(row[1], []).append([row[0], row[2], first_time[row[1]], row[7]])

    return len(tmp_list), d

def get_sample_project(d):
    return d.values()[0][0][3]  # the first one, for instance

