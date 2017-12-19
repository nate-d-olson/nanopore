import collections
import glob
import os
import pprint
import sys

from admiral import jobmanagers
from admiral import remote
from admiral import slurm

import basecalling
import experiments
import fast5_archives
import mapping

# import slurm

BASE_PATH = "/oak/stanford/groups/msalit/nspies/nanopore"

def run_dir(run_name):
    return f"{BASE_PATH}/raw/{run_name}"

def fast5_dir(run_name):
    return f"{run_dir(run_name)}/fast5"

def fastq_dir(run_name):
    return f"{run_dir(run_name)}/fastq"

def mappings_dir(run_name):
    return f"{run_dir(run_name)}/aln"

def combined_bam(run_name):
    return f"{mappings_dir(run_name)}/{run_name}.combined.sorted.bam"


def _jobmanager():
    return slurm.SLURM_Jobmanager(batch_dir="output", log_dir="output")

def iter_runs(metadata):
    for flowcell_name, flowcell_info in metadata.items():
        for runinfo in flowcell_info["datasets"]:
            run_name = runinfo["name"]
            
            yield flowcell_name, flowcell_info, runinfo, run_name    

def basecalling_complete(metadata):
    for _, outpath, _, _ in get_basecalling_args(metadata, 0):
        if not os.path.exists(outpath):
            return False
    return True

def mapping_complete(metadata):
    for _,_,_,run_name in iter_runs(metadata):
        if not os.path.exists(combined_bam(run_name)):
            return False
    return True

### Archiving

def archive_run(run_base_path, run_name):
    chunks = glob.glob(f"{run_base_path}/*")
    chunks = [chunk for chunk in chunks if os.path.isdir(chunk)]

    for chunk in chunks:
        chunk_name = os.path.basename(chunk)
        tar_path = f"{run_base_path}/{run_name}_{chunk_name}.tar"

        yield chunk, tar_path

def get_archiving_args(metadata):
    args = []
    for _,_,run_info,run_name in iter_runs(metadata):
        for chunk, tar_path in archive_run(fast5_dir(run_name), run_name):
            should_remove = run_info["completed"]
            args.append([chunk, tar_path, should_remove])
    return args

def launch_archiving(metadata):
    args = get_archiving_args(metadata)
    print(args)
    njobs = min(len(args), 128)
    job = remote.run_remote(
        fast5_archives.archive_chunk, _jobmanager(),
        job_name="archive_fast5s", args=args, job_dir="output",
        overwrite=True, njobs=njobs, queue="owners", mem="8g")

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))


### Baseacalling

def get_basecalling_args(metadata, threads):
    args = []
    for _,flowcell_info,runinfo,run_name in iter_runs(metadata):      
        os.makedirs(fastq_dir(run_name), exist_ok=True)

        for fast5_archive in glob.glob(f"{fast5_dir(run_name)}/{run_name}_*.tar"):
            chunk_name = os.path.splitext(os.path.basename(fast5_archive))[0]+".fastq.gz"
            outpath = f"{fastq_dir(run_name)}/{chunk_name}"

            config = {"flowcell":flowcell_info["flowcell_type"],
                      "kit":runinfo["kit"]}

            args.append([fast5_archive, outpath, config, threads]) 

    return args

def launch_basecalling(metadata):
    threads = 8
    args = get_basecalling_args(metadata, threads)
    if len(args) == 0:
        print("No directories found for basecalling...")
        return
    
    njobs = min(len(args), 128)
    if njobs < len(args):
        # should use int(math.ceil(len(args)/128))) or something like that
        raise Exception("NEED TO ADJUST RUN TIME accordingly")

    job = remote.run_remote(
        basecalling.run_basecalling_locally, _jobmanager(),
        job_name="basecalling", args=args, job_dir="output",
        overwrite=True, njobs=njobs, queue="owners", 
        cpus=threads, mem=f"{8*threads}g")

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))


### Mapping

def get_mapping_args(metadata, threads):
    args = []
    for _,_,_,run_name in iter_runs(metadata):
        os.makedirs(mappings_dir(run_name), exist_ok=True)

        for fastq in glob.glob(f"{fastq_dir(run_name)}/*.fastq.gz"):
            chunk_name = os.path.splitext(os.path.basename(fastq))[0]+".sorted.bam"
            outpath = f"{mappings_dir(run_name)}/{chunk_name}"
            args.append([fastq, outpath, threads])

    return args

def _get_merge_bams_args(run_name):
    bams = []
    for bam in glob.glob(f"{mappings_dir(run_name)}/*.sorted.bam"):
        if "combined" in bam: continue

        bams.append(bam)

    return [combined_bam(run_name), bams]

def get_merge_bams_args(metadata):
    args = []
    for _,_,_,run_name in iter_runs(metadata):
        cur_args = _get_merge_bams_args(run_name)
        if len(cur_args[1]) == 0:
            print(f"No bams to merge {run_name}")
            continue
        args.append(cur_args)

    return args

def launch_mapping(metadata):
    threads = 4
    args = get_mapping_args(metadata, threads)

    if len(args) == 0:
        print("No fastq files found for mapping...")
        return
        
    njobs = min(len(args), 128)
    if njobs < len(args):
        # should use int(math.ceil(len(args)/128))) or something like that
        raise Exception("NEED TO ADJUST RUN TIME accordingly")

    job = remote.run_remote(
        mapping.run_mapping, _jobmanager(),
        job_name="mapping", args=args, job_dir="output",
        overwrite=True, njobs=njobs, queue="owners", mem="48g", cpus=threads)

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))


def launch_merge_bams(metadata):
    args = get_merge_bams_args(metadata)

    if len(args) == 0:
        print("No bam files found to merge...")
        return
        
    njobs = min(len(args), 128)
    if njobs < len(args):
        # should use int(math.ceil(len(args)/128))) or something like that
        raise Exception("NEED TO ADJUST RUN TIME accordingly")

    job = remote.run_remote(
        mapping.merge_bams, _jobmanager(),
        job_name="merge_bams", args=args, job_dir="output",
        overwrite=True, njobs=njobs, queue="owners", mem="8g")

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))



def filter_completed_datasets(metadata):
    to_run = collections.OrderedDict()
    completed = []

    for name, info in metadata.items():
        if basecalling_complete({name:info}) and mapping_complete({name:info}):
            completed.append(name)
        else:
            to_run[name] = info

    return to_run, completed
    
def main():      
    metadata = experiments.load_experiment_metadata()
    if len(sys.argv) == 2:
        flow_cell_id = sys.argv[1]
        metadata = {flow_cell_id: metadata[flow_cell_id]}
    else:
        metadata, completed = filter_completed_datasets(metadata)

        if len(completed) > 0:
            print("-"*50)
            print(f"Already completed: {','.join(completed)}")
            print("...skipping.")
            print("-"*50)

    pprint.pprint(metadata)


    launch_archiving(metadata)
    launch_basecalling(metadata)
    launch_mapping(metadata)
    launch_merge_bams(metadata)


if __name__ == '__main__':
    main()
