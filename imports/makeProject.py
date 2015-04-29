# Create temporary directory for each project
def makeProject(inputPath):
    import os
    from os.path import normpath, basename
    dirName = basename(normpath(inputPath))
    os.makedirs(dirName)
    os.makedirs(dirName+"/pbs")

def getDirName(inputPath):
    import os
    from os.path import normpath, basename
    dirName = basename(normpath(inputPath))
    return dirName


