import subprocess


if __name__ == '__main__':
    # submit the job to the queue
    p = subprocess.Popen(
        ['qsub', 'pbs/run.pbs'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    # gather output and errors
    out, err = p.communicate()
    print(out)
