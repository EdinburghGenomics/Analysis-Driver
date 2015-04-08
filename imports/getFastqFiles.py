# Return a list with the names of fastq files
# within the given directory
def getFastqFiles(path):
    import os
    # TODO: Check the paths exists
    listNames=[]
    for file in os.listdir(path):
        if file.endswith(".fastq"):
            listNames.append(file)

    return listNames
