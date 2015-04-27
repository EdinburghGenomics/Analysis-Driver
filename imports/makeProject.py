# Create temporary directory for each project
def makeProject(inputPath):
    from os.path import normpath, basename
    dirName = basename(normpath(inputPath))
    os.makedirs(dirName)
    
    os.makedirs(dirName+"/pbs")
