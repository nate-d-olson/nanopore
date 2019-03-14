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

# BASE_PATH = "/scratch/groups/msalit/nanopore"
## Path for testing
BASE_PATH = "/scratch/groups/msalit/nanopore/nanopore-test"

def run_dir(run_name):
    return f"{BASE_PATH}/raw/{run_name}"

def fast5_dir(run_name):
    return f"{run_dir(run_name)}/fast5"

def fastq_dir(run_name):
    return f"{run_dir(run_name)}/fastq"

def mappings_dir(run_name, ref_name):
    return f"{run_dir(run_name)}/aln_{ref_name}"

def combined_bam(run_name, ref_name):
    return f"{mappings_dir(run_name, ref_name)}/{run_name}.combined.sorted.bam"


def _jobmanager():
    return slurm.SLURM_Jobmanager(batch_dir="output", log_dir="output")

def iter_runs(metadata, exclude_basecalled=False):
    for flowcell_name, flowcell_info in metadata.items():
        for runinfo in flowcell_info["datasets"]:
            run_name = runinfo["name"]
            
            yield flowcell_name, flowcell_info, runinfo, run_name    

def basecalling_complete(metadata):
    for _, outpath, _, _ in get_basecalling_args(metadata, 0):
        if not os.path.exists(outpath):
            print("Not found:", outpath)
            return False
    return True

def mapping_complete(metadata, ref_name):
    for _,_,_,run_name in iter_runs(metadata):
        if not os.path.exists(combined_bam(run_name, ref_name)):
            print("Not found:", combined_bam(run_name, ref_name))
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
    for _,_,run_info,run_name in iter_runs(metadata, exclude_basecalled=True):
        for chunk, tar_path in archive_run(fast5_dir(run_name), run_name):
            should_remove = run_info["completed"]
            args.append([chunk, tar_path, should_remove])
    return args

def launch_archiving(metadata):
    args = get_archiving_args(metadata)

    if len(args) == 0:
        print("Archiving complete; skipping...")
        return
    
    njobs = min(len(args), 128)
    job = remote.run_remote(
        fast5_archives.archive_chunk, _jobmanager(),
        job_name="archive_fast5s", args=args, job_dir="output",
        overwrite=True, njobs=njobs, queue="owners", mem="8g")

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))


### Basecalling
def get_guppy_config(flowcell, kit, platform):
    """
    defines config file for use in basecalling based on flowcell type, sequencing kit, and platform
    """

    ## Promethion
    if platform == "promethion":
        return "dna_r9.4.1_450bps_flipflop_prom.cfg"
    elif  flowcell == "FLO-MIN107" and kit == "SQK-RAD003":
        return "dna_r9.5_450bps.cfg"
    elif  flowcell == "FLO-MIN106" and kit == "SQK-RAD003":
        return "dna_r9.4.1_450bps_flipflop.cfg"
    elif  flowcell == "FLO-MIN106" and kit == "SQK-RAD004":
        return "dna_r9.4.1_450bps_flipflop.cfg"
    else:
        print(f"Config file not defined for {kit}, {flowcell}, and {platform} combination.")
        return "No defined config"

def get_basecalling_args(metadata, threads):
    args = []
    for _,flowcell_info,runinfo,run_name in iter_runs(metadata, exclude_basecalled=True):
        os.makedirs(fastq_dir(run_name), exist_ok=True)

        ## List of tar'd fast5s or multi-seq fast5s
        fast5s = glob.glob(f"{fast5_dir(run_name)}/*{run_name}_*.tar")

        for fast5_archive in fast5s:
            chunk_name = os.path.splitext(os.path.basename(fast5_archive))[0]+".fastq.gz"
            outpath = f"{fastq_dir(run_name)}/{chunk_name}"

            guppy_config = get_guppy_config(flowcell = flowcell_info["flowcell_type"], 
                                            kit = runinfo["kit"], 
                                            platform = runinfo["platform"])

            config = {"flowcell":flowcell_info["flowcell_type"],
                      "kit":runinfo["kit"],
                      "config":guppy_config}

            if not os.path.exists(outpath) or (os.path.getmtime(outpath) < os.path.getmtime(fast5_archive)):
                args.append([fast5_archive, outpath, config, threads]) 

    return args

def launch_basecalling(metadata):
    threads = 16
    args = get_basecalling_args(metadata, threads)
    if len(args) == 0:
        print("No directories found for basecalling...")
        return

    print(f"Running basecalling on {len(args)} fast5 files...")
    
    njobs = min(len(args), 500)
    
    if njobs > len(args):
        print("MESSAGE: Number of tar balls > 500. Only processing first 500. Rerun pipeline to process next batch of tars.")

    job = remote.run_remote(
        basecalling.run_basecalling_locally, _jobmanager(),
        job_name="basecalling", args=args[0:njobs], job_dir="output",
        overwrite=True, njobs=njobs, queue="msalit,owners,normal", 
        cpus=threads, mem=f"{8+threads}g", time="4h")

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))


### Mapping
def get_genome_path(ref_name):
    if ref_name == "hs37d5":
        return "/oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa"
    elif ref_name == "GRCh38":
        return "/oak/stanford/groups/msalit/shared/genomes/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
    else:
        print("Invalid ref_name")
        return ""

def get_mapping_args(metadata, ref_name, threads):
    args = []
    
    genome_path = get_genome_path(ref_name)
    
    for flowcell_id,flowcell_info,runinfo,run_name in iter_runs(metadata):
        os.makedirs(mappings_dir(run_name, ref_name), exist_ok=True)
        
        ## Read Group Args
        guppy_config = get_guppy_config(flowcell = flowcell_info["flowcell_type"], 
                                        kit = runinfo["kit"], 
                                        platform = runinfo["platform"])

        read_group_args = {"run":run_name,
                            "flowcell_id":flowcell_id,
                            "platform_model":runinfo["platform"],
                            "sample":runinfo["sample"],
                            "date":runinfo["date"],
                            "guppy_config":guppy_config,
                            "flowcell_type":flowcell_info["flowcell_type"],
                            "flowcell_kit":runinfo["kit"]}

        for fastq in glob.glob(f"{fastq_dir(run_name)}/*.fastq.gz"):
            chunk_name = os.path.splitext(os.path.basename(fastq))[0]+".sorted.bam"
            outpath = f"{mappings_dir(run_name, ref_name)}/{chunk_name}"
            args.append([fastq, outpath, read_group_args, genome_path, threads])

    return args

def _get_merge_bams_args(run_name, ref_name):
    bams = []
    for bam in glob.glob(f"{mappings_dir(run_name, ref_name)}/*.sorted.bam"):
        if bam.endswith("combined.sorted.bam"): continue

        bams.append(bam)

    return [combined_bam(run_name, ref_name), bams]

def get_merge_bams_args(metadata, ref_name):
    args = []
    for _,_,_,run_name in iter_runs(metadata):
        cur_args = _get_merge_bams_args(run_name, ref_name)
        if len(cur_args[1]) == 0:
            print(f"No bams to merge {run_name}")
            continue
        args.append(cur_args)

    return args

def launch_mapping(metadata, ref_name):
    threads = 4
    args = get_mapping_args(metadata, ref_name, threads)
    print(f"Aligning {len(args)} fastq files...")
    
    if len(args) == 0:
        print("No fastq files found for mapping...")
        return
        
    njobs = min(len(args), 500)

    job = remote.run_remote(
        mapping.run_mapping, _jobmanager(),
        job_name="mapping", args=args[0:njobs], job_dir="output",
        overwrite=True, njobs=njobs, queue="owners,msalit,normal", mem="48g", cpus=threads)

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))


def launch_merge_bams(metadata, ref_name):
    args = get_merge_bams_args(metadata, ref_name)

    if len(args) == 0:
        print("No bam files found to merge...")
        return
        
    njobs = min(len(args), 128)

    job = remote.run_remote(
        mapping.merge_bams, _jobmanager(),
        job_name="merge_bams", args=args, job_dir="output",
        overwrite=True, njobs=njobs, queue="owners,msalit,normal", cpus=4, mem="24g")

    print(jobmanagers.wait_for_jobs([job], progress=True, wait=5.0))



def filter_completed_datasets(metadata):
    to_run = collections.OrderedDict()
    completed = []

    for name, info in metadata.items():
        if basecalling_complete({name:info}) and mapping_complete({name:info}, "hs37d5") and mapping_complete({name:info}, "GRCh38"):
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

    if not basecalling.is_guppy_installed():
        print("[ERROR] ONT basecallier Guppy) must \n"
              "[ERROR] be installed and accessible from the current environment.")
        sys.exit(1)

    print("Archiving...")
    launch_archiving(metadata)

    print("Basecalling...")
    launch_basecalling(metadata)

    print("Mapping hs37d5...")
    launch_mapping(metadata, ref_name = "hs37d5")

    print("Mapping GRCh38...")
    launch_mapping(metadata, ref_name = "GRCh38")

#    print("Merging bams...")
#    launch_merge_bams(metadata)


if __name__ == '__main__':
    main()
