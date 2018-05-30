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

def pull_raw_data(fast5_archive_path):
    """
    untars a fast5 archive onto the local hard-drive
    returns the path to the decompressed directory
    """

    local_dir = os.environ["LOCAL_SCRATCH"]
    temp_dir = tempfile.TemporaryDirectory(dir=local_dir)
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
    
    command = "read_fast5_basecaller.py \
                      --flowcell {flowcell} \
                      --kit {kit} \
                      -o fastq -q 0 \
                      -t {threads} \
                      -s {outdir} \
                      --basecaller.max_events=10000" # this is the default in the latest version of albacore

    command = command.format(
        flowcell=config["flowcell"],
        kit=config["kit"],
        threads=threads,
        outdir=workingdir.name)

    stdin = "\n".join(fast5_paths)
    
    print(command)
    
    process = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
    process.communicate(input=stdin.encode())
    
    pass_fastqs = glob.glob(f"{workingdir.name}/workspace/pass/*.fastq")
    fail_fastqs = glob.glob(f"{workingdir.name}/workspace/fail/*.fastq")

    #assert len(pass_fastqs) <= 1
    #assert len(fail_fastqs) <= 1

    pass_count = 0
    if len(pass_fastqs):
        pass_count = wc(pass_fastqs[0])

    fail_count = 0
    if len(fail_fastqs):
        fail_count = wc(fail_fastqs[0])

    print(f"Pass: {pass_count}   Fail: {fail_count}")

    if len(pass_fastqs) == 0:
        open(out_fastq_gz, "w")
        return
    
    assert process.returncode == 0
    # shutil.move(pass_fastqs[0], out_fastq)


    subprocess.check_call(f"pigz -c {' '.join(pass_fastqs)} > {out_fastq_gz}", shell=True)

    shutil.move(f"{workingdir.name}/sequencing_summary.txt", f"{out_fastq_gz}.sequencing_summary.txt") # to test

def run_basecalling_locally(fast5_archive_path, out_fastq, config, threads=1):
    print("extracting data to local storage...")
    fast5_dir = pull_raw_data(fast5_archive_path)

    fast5_paths = glob.glob(f"{fast5_dir.name}/*/*.fast5")#[:50]
    if len(fast5_paths) == 0:
        fast5_paths = glob.glob(f"{fast5_dir.name}/*.fast5")#[:50]
        
    print(f"performing basecalling on {len(fast5_paths)} reads...")
    perform_basecalling(fast5_paths, out_fastq, config, threads)



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
