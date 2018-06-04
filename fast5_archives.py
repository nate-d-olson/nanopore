# import glob
import os
import shutil
import subprocess
import sys
import tarfile
import time

class CombineTars:
    def __init__(self, outpath):
        self.out = tarfile.TarFile(outpath, "w")
        self.subdirs = set()

    def append_contents(self, other_tar_path, dir_name):
        print(f"  {other_tar_path} - {dir_name}")
        other_tar = tarfile.TarFile(other_tar_path)

        subdir = tarfile.TarInfo(f"./{dir_name}")
        if subdir.name not in self.subdirs:
            subdir.type = tarfile.DIRTYPE
            subdir.mode = 0o755
            subdir.mtime = time.time()
            self.out.addfile(subdir)
            self.subdirs.add(subdir.name)
            
        for member in other_tar:
            if member.name.endswith(".fast5"):
                member_name = member.name.split("/")[-1]
                dest_info = tarfile.TarInfo(f"./{dir_name}/{member_name}")
                dest_info.mtime = member.mtime
                dest_info.size = member.size
                dest_info.mode = member.mode
                dest_info.type = member.type
                self.out.addfile(dest_info, other_tar.extractfile(member))


def archive_chunk(fast5_dir, tar_path, remove=False):
    """
    uses tar to archive a directory of fast5 files
    """
    top_dir, _, chunk_dir = fast5_dir.rstrip("/").rpartition("/")

    command = f"tar -cf {tar_path} -C {top_dir} {chunk_dir}"
    print(f"Running '{command}'...")
    subprocess.check_call(command, shell=True)

    if not validate_tar(fast5_dir, tar_path):
        raise Exception(f"Failed to completely tar {fast5_dir} into {tar_path}")

    if remove:
        shutil.rmtree(fast5_dir)
        

def validate_tar(original_dir, tar_path):
    tar = tarfile.TarFile(tar_path)
    tar_contents = [os.path.basename(x) for x in tar.getnames()]

    original_files = os.listdir(original_dir)

    return tar_contents <= original_files


# def archive_run(run_base_path, run_name):
#     print("*"*30)
#     chunks = glob.glob(f"{run_base_path}/*")
#     chunks = [chunk for chunk in chunks if os.path.isdir(chunk)]

#     for chunk in chunks:
#         chunk_name = os.path.basename(chunk)
#         tar_path = f"{run_base_path}/{run_name}_{chunk_name}.tar"

#         archive_chunk(chunk, tar_path)

# def test():
#     archive_run(
#         "/oak/stanford/groups/msalit/nspies/nanopore/raw/20170921_2217_HG002_ultralong_plug_b/fast5/",
#         "20170921_2217_HG002_ultralong_plug_b")

if __name__ == '__main__':
    # archive_run(sys.argv[1], sys.argv[2])
    ct = CombineTars(sys.argv[1])
    for other_tar in sys.argv[2:]:
        path, name = other_tar.split(",")
        ct.append_contents(path, name)
