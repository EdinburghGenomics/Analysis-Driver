# Create temporary directory for each project
def makeProject(inputPath):
    import os
    os.makedirs(inputPath)
    os.makedirs(inputPath+"/pbs")



def getDirName(inputPath):
    import os
    from os.path import normpath, basename
    dirName = basename(normpath(inputPath))
    return str(dirName)
