#!/usr/bin/bash
#SBATCH -n 8 # number of cores
#SBATCH --mem 48g # memory pool for all cores
#SBATCH --job-name=cut_release37
#SBATCH --time=1-00:00:00
#SBATCH --partition=msalit
#SBATCH --mail-type=ALL

export PATH=$PATH:/oak/stanford/groups/msalit/ndolson/nanopore-pipeline/

## Loading required modules
module load python/3.6.1 R biology samtools

cd /oak/stanford/groups/msalit/ndolson

## activating env
source nanopore_env/bin/activate

## running run release script
bash nanopore-pipeline/cut_release.sh ultra-long-ont

## deactivating env
deactivate
