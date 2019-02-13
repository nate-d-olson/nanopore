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

def is_albacore_installed():
    try:
        print("Checking albacore (ONT basecaller) is installed and can be run...")
        subprocess.check_call("read_fast5_basecaller.py -v", shell=True)
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

    # if we need the subdirectory name, probably unnecessary
    # base = fast5_archive.members[0].name

    fast5_archive.extractall(path=temp_dir.name)
        
    return temp_dir



def perform_basecalling(fast5_paths, out_fastq_gz, config, threads):
    """
    basecalls the specified fast5 files into out_fastq
    """

    workingdir = tempfile.TemporaryDirectory()

    ## Albacore command
    # command = "read_fast5_basecaller.py \
    #                   --flowcell {flowcell} \
    #                   --kit {kit} \
    #                   -o fastq -q 0 \
    #                   -t {threads} \
    #                   -s {outdir} \
    #                   --basecaller.max_events=10000" # this is the default in the latest version of albacore

    ## Guppy command V2.3.1
    command = "guppy_basecaller \
		            --input_path {fast5_path} \
                    --save_path {outdir} \
                    --flowcell {flowcell} \
                    --kit {kit} \
                    --num_callers {threads}"

    command = command.format(
        fast5_path=fast5_paths,
        outdir=workingdir.name,
        flowcell=config["flowcell"],
        kit=config["kit"],
        threads=threads)

    ## Guppy command V2.3.1 with flip-flop for flowcell: FLO-MIN106 and kit:SQK-RAD004

    command = "guppy_basecaller \
		            --input_path {fast5_path} \
                    --save_path {outdir} \
                    -c dna_r9.4.1_450bps_flipflop.cfg \
                    --num_callers {threads}"

    command = command.format(
        fast5_path=fast5_paths,
        outdir=workingdir.name,
        threads=threads)

    ## Albacore
    # stdin = "\n".join(fast5_paths)

    
    # process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
    # process.communicate(input=stdin.encode())

    ## For guppy not using standard input to pass list of fast5s
    print("Guppy command")
    print(command)
    subprocess.run(command, shell = True)

    # pass_fastqs = glob.glob(f"{workingdir.name}/workspace/pass/*.fastq")
    out_fastqs = glob.glob(f"{workingdir.name}/*.fastq")
    fail_fastqs = glob.glob(f"{workingdir.name}/workspace/fail/*.fastq")

    # pass_count = 0
    # if len(pass_fastqs):
    #     pass_count = wc(pass_fastqs[0])

    # fail_count = 0
    # if len(fail_fastqs):
    #     fail_count = wc(fail_fastqs[0])

    # print(f"Pass: {pass_count}   Fail: {fail_count}")
    # print(f"Pass: {pass_count}   Fail: {fail_count}")

    # # if len(pass_fastqs) == 0:
    # #     open(out_fastq_gz, "w")
    # #     return
    
    # # assert process.returncode == 0
    print("Output fastqs")
    print(out_fastqs)
    # subprocess.check_call(f"pigz -c {' '.join(pass_fastqs)} > {out_fastq_gz}", shell=True)
    subprocess.check_call(f"pigz -c {' '.join(out_fastqs)} > {out_fastq_gz}", shell=True)

    shutil.move(f"{workingdir.name}/sequencing_summary.txt", f"{out_fastq_gz}.sequencing_summary.txt") # to test


def run_basecalling_locally(fast5_archive_path, out_fastq, config, threads):
    print("extracting data to local storage...")
    fast5_dir = pull_raw_data(fast5_archive_path)

    ## Albacore - fast5 files
    fast5_paths = glob.glob(f"{fast5_dir.name}/*.fast5")
 
    print(f"performing basecalling on {len(fast5_paths)} reads...")
    perform_basecalling(fast5_dir.name, out_fastq, config, threads)



# def test():
#     print("extracting data to local storage")
#     fast5_archive = "/oak/stanford/groups/msalit/nspies/nanopore/raw/20170927_1910_HG002_ultralong_plug_filter/fast5/20170927_1910_HG002_ultralong_plug_filter_0.tar"
#     fast5_dir = pull_raw_data(fast5_archive)

#     print(fast5_dir)

#     fast5_paths = glob.glob(f"{fast5_dir}/*/*.fast5")[:10]
#     out_fastq = "asdf.fastq"

#     configuration = {"flowcell":"FLO-MIN106",
#                      "kit":"SQK-RAD003",
#                      "threads":1}

#     print(f"performing basecalling on {len(fast5_paths)} reads...")
#     perform_basecalling(fast5_paths, out_fastq, configuration)


if __name__ == '__main__':
    import json, sys
    
    config = json.loads(sys.argv[3])
    print(config)
    
    run_basecalling_locally(sys.argv[1], sys.argv[2], config, config.get("threads"))
