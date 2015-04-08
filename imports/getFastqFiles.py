# Return a list with the names of fastq files
# within the given directory
def getFastqFiles(path):
    import os,sys
    # TODO: Check the paths exists
    try:
        assert ( os.path.isdir(path) == True)
    except AssertionError:
        print "Expected a valid directory to find fastq files"
        sys.exit()
    listNames=[]
    for file in os.listdir(path):
        if file.endswith(".fastq"):
            listNames.append(file)
    return listNames
