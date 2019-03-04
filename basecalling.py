"""
Untars fast5 files to local storage, performs basecalling with albacore, then puts 
the resulting fastq and sequencing_summary files at the specified location in
shared storage.
"""

import glob
import tempfile
import os
import shutil
import subprocess
import tarfile


def wc(path):
    count = 0
    for line in open(path):
        count += 1
    return count

def is_guppy_installed():
    try:
        print("Checking guppy (ONT basecaller) is installed and can be run...")
        subprocess.check_call("guppy_basecaller -v", shell=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False

def pull_raw_data(fast5_archive_path):
    """
    untars a fast5 archive onto the local hard-drive
    returns the path to the decompressed directory
    """

    local_dir = os.environ["LOCAL_SCRATCH"]
    temp_dir = tempfile.TemporaryDirectory(dir=local_dir)

    ## Check if archived fast5s or directory of fast5s
    fast5_archive = tarfile.TarFile(fast5_archive_path)

    fast5_archive.extractall(path=temp_dir.name)
        
    return temp_dir



def perform_basecalling(fast5_paths, out_fastq_gz, config, threads):
    """
    basecalls the specified fast5 files into out_fastq
    """

    workingdir = tempfile.TemporaryDirectory()

    ## Guppy command V2.3.1 with flip-flop for flowcell: FLO-MIN106 and kit:SQK-RAD004
    command = "guppy_basecaller \
		        --input_path {fast5_path} \
		        --recursive \
                --save_path {outdir} \
                -c {config} \
                --num_callers {threads} \
		        --cpu_threads_per_caller 4"

    command = command.format(
        fast5_path=fast5_paths,
        outdir=workingdir.name,
        config=config["config"],
        threads=threads//4)

    ## For guppy not using standard input to pass list of fast5s
    subprocess.run(command, shell = True)

    out_fastqs = glob.glob(f"{workingdir.name}/*.fastq")
    
    subprocess.check_call(f"pigz -c {' '.join(out_fastqs)} > {out_fastq_gz}", shell=True)

    ## Seq summary 
    shutil.move(f"{workingdir.name}/sequencing_summary.txt", f"{out_fastq_gz}.sequencing_summary.txt") 
    
    ## Basecalling log file
    guppy_log=glob.glob(f"{workingdir.name}/*.log")
    if len(guppy_log) == 1:
        shutil.move(guppy_log[0], f"{out_fastq_gz}.guppy_basecaller_log.log") 
    elif len(guppy_log) > 1:
        print("More than one guppy basecalling log file - concatenating")
        print(guppy_log)
        with open(f"{out_fastq_gz}.guppy_basecaller_log.log", 'w') as outfile:
            for fname in guppy_log:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    else:
        print("no basecalling log files found")

    ## Guppy telemetry file
    shutil.move(f"{workingdir.name}/sequencing_telemetry.js", f"{out_fastq_gz}.sequencing_telemetry.js") 

def run_basecalling_locally(fast5_archive_path, out_fastq, config, threads):
    print("extracting data to local storage...")
    fast5_dir = pull_raw_data(fast5_archive_path)

    ## Albacore - fast5 files
    fast5_paths = glob.glob(f"{fast5_dir.name}/*.fast5")
 
    print(f"performing basecalling on {len(fast5_paths)} reads...")
    perform_basecalling(fast5_dir.name, out_fastq, config, threads)


if __name__ == '__main__':
    import json, sys
    
    config = json.loads(sys.argv[3])
    print(config)
    
    run_basecalling_locally(sys.argv[1], sys.argv[2], config, config.get("threads"))
