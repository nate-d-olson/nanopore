#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 18g # memory pool for all cores
#SBATCH --job-name=rel_fqss
#SBATCH --time=1-00:00:00
#SBATCH --partition=msalit
#SBATCH --mail-type=ALL


## Loading required modules
module load python/3.6.1 R biology samtools

cd /oak/stanford/groups/msalit/ndolson

## activating env
source nanopore_env/bin/activate

## Adding pipeline code dir to path
export PATH=$PATH:/oak/stanford/groups/msalit/ndolson/nanopore-pipeline/

## running run release script
bash nanopore-pipeline/cut_release.sh ultra-long-ont

## deactivating env
deactivate
