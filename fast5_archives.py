import glob
import os
import subprocess
import sys

def archive_chunk(fast5_dir, tar_path):
    """
    uses tar to archive a directory of fast5 files
    """
    top_dir, _, chunk_dir = fast5_dir.rstrip("/").rpartition("/")

    command = f"tar -cf {tar_path} -C {top_dir} {chunk_dir}"
    print(f"Running '{command}'...")
    subprocess.check_call(command, shell=True)


# def archive_run(run_base_path, run_name):
#     print("*"*30)
#     chunks = glob.glob(f"{run_base_path}/*")
#     chunks = [chunk for chunk in chunks if os.path.isdir(chunk)]

#     for chunk in chunks:
#         chunk_name = os.path.basename(chunk)
#         tar_path = f"{run_base_path}/{run_name}_{chunk_name}.tar"

#         archive_chunk(chunk, tar_path)

def test():
    archive_run(
        "/oak/stanford/groups/msalit/nspies/nanopore/raw/20170921_2217_HG002_ultralong_plug_b/fast5/",
        "20170921_2217_HG002_ultralong_plug_b")

if __name__ == '__main__':
    archive_run(sys.argv[1], sys.argv[2])