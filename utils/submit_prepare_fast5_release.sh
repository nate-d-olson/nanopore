#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 18g # memory pool for all cores
#SBATCH --job-name=rel_fast5
#SBATCH --time=2-00:00:00
#SBATCH --partition=msalit
#SBATCH --mail-type=ALL


## Loading required modules
module load python/3.6.1 R biology samtools

cd $OAK/ndolson

## activating env
source nanopore_env/bin/activate

## Adding pipeline code dir to path
export PATH=$PATH:$OAK/ndolson/nanopore-pipeline/

## Moving to processing directory
cd nanopore-pipeline

## running run release script
python3 prepare_fast5_release.py

## deactivating env
deactivate
