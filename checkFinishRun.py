mypath = "/home/U008/lcebaman/scripts/data/Finish.txt"

#fileExists(mypath)


#def fileExists(file)
import os.path,sys
try:
    print os.path.isfile(mypath) 
except os.error:
    print 'ERROR in fileExists'

sys.exit()
