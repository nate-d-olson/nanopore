import subprocess

def launch_job(command, time=200, cores=1, mem=None, queue="msalit"):
    if not command.startswith("#!"):
        command = "#!/bin/bash\n" + command

    if mem is None:
        mem = cores * 8
    assert mem < 256, "mem should be in GB"
    
    print("Launching: '{}'".format(command))
    # time in minutes
    batch_command = f"echo '{command}' | sbatch -t {time} -N 1 -c {cores} --mem {mem*1000} -p {queue}"

    batch_proc = subprocess.run(batch_command, shell=True, stdout=subprocess.PIPE)

    stdout = batch_proc.stdout.decode("utf-8")

    print(f"sbatch result: '{stdout}'")
    assert stdout.startswith("Submitted")
    jobid = stdout.splitlines()[0].split()[-1]
    return jobid

def test():
    print("::", launch_job("echo 12345"))

if __name__ == '__main__':
    test()