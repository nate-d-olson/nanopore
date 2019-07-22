
# coding: utf-8

# # Fast5 tar archives for dataset release

# Next Steps:
#  - Testing output
#  - work into script for job submission

# In[62]:


import os
import pprint
import sys
import collections
import glob

import experiments

## Moving to raw data directory
os.chdir("/scratch/groups/msalit/nanopore/raw")

from fast5_archives import CombineTars, natural_sort


# ## Plan
# 1. Use experiments.py to load metadata
# 2. Generate list of directories for each flowcell
# 3. Pass directory list to a bash script with the combine_fast5s code, Script takes flowcell and a list of directories as input, Outputs tar archives in 50 Gb chunks!

# ## Loading Metadata

# In[3]:


metadata = experiments.load_experiment_metadata()


# ## Generating list of directories for each flowcell

# In[34]:


## Getting list for directories by flowcell
def run_dir(run_name):
    return f"{run_name}"

def fast5_dir(run_name):
    return f"{run_dir(run_name)}/fast5"

def get_flowcell_runs(metadata):
    flowcell_dict = collections.defaultdict(dict)
    
    for flowcell_id, flowcell_info in metadata.items():
        flowcell_dict[flowcell_id] = []
        for runinfo in flowcell_info["datasets"]:
            flowcell_dict[flowcell_id].append(runinfo["name"])
    
    return(flowcell_dict)

flowcell_dict = get_flowcell_runs(metadata)


# ## Archiving fast5s

# In[73]:


release_path = f"../release/ultra-long-ont-fast5s"
def get_tree_size(path):
    """Return total size of files in given path and subdirs."""
    total = 0
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            total += get_tree_size(entry.path)
        else:
            total += entry.stat(follow_symlinks=False).st_size
    return total

def combine_tars(tar_filename, tar_list):
    ct = CombineTars(tar_filename)
    for other_tar in tar_list:
        path, name = other_tar.split(",")
        ct.append_contents(path, name)
    
n = 1  
total_flowcells = len(flowcell_dict)
for flowcell, flowcell_runs in flowcell_dict.items():
    print("Processing {}: flowcell {} of {}".format(flowcell, n, total_flowcells))
    
    ## Fast5 index number
    idx_fast5_list = []
    idx_total_size = 0
    idx = 0
    
    ## Iterating through flowcell runs
    for fc_run in flowcell_runs:
        fc_dir = fast5_dir(fc_run)
        
        for fast5 in glob.glob(f"{fc_dir}/*{fc_run}_*.tar"):
            fast5_size = os.path.getsize(fast5)
            if fast5_size == 0:
                print(f"Empty fast5 tar: {fast5}")
            else:
                if idx_total_size + fast5_size > 50*1e9:
                    combine_tars(f"{release_path}/{flowcell}_{idx}.tar", fast5_list)
                    ## Resetting total size and fast5 list 
                    idx += 1
                    idx_total_size = fast5_size
                    idx_fast5_list = [f"{fast5},{fc_run}"]
                else:
                    idx_fast5_list.append(f"{fast5},{fc_run}")
                    idx_total_size += fast5_size
        
        ## Archiving Remaining fast5s
        combine_tars(f"{release_path}/{flowcell}_raw_fast5s_{idx}.tar", idx_fast5_list)
    n += 1      

