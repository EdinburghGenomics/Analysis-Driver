"""
To run shell commands
"""
import subprocess
                                             
# submit the job to the queue                
p = subprocess.Popen(['qsub', 'pbs/run.pbs'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# gather output and errors                     
out, err = p.communicate()
print(out)
