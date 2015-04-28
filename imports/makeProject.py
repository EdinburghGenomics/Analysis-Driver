# Create temporary directory for each project
def makeProject(inputPath):
    import os
    from os.path import normpath, basename
    dirName = basename(normpath(inputPath))
    os.makedirs(dirName+"_out")
    os.makedirs(dirName+"_out"+"/pbs")
    return dirName


