# Check that if file has been created 
# Input: Path to the filename + name of the file
# Output: True(file exists) or False (file has not been created yet)
def fileExists(file):
    import os.path,sys
    try:
        exists = os.path.isfile(file) 
    except os.error:
        print 'ERROR in fileExists'
        sys.exit()
   
    return exists

# Unit Test
filename = "/home/U008/lcebaman/scripts/data/Finish.txt"
print fileExists(filename)
