#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --job-name=ONT-test
#SBATCH --time=10:00:00
#SBATCH --partition=msalit,owners
#SBATCH --mail-type=ALL

BASEPATH=/oak/stanford/groups

## Loading required modules
module load python/3.6.1 R biology samtools

export PATH=$PATH:${BASEPATH}/msalit/ndolson/minimap2/
export PATH=$PATH:${BASEPATH}/msalit/ndolson/pigz/
export PATH=$PATH:${BASEPATH}/msalit/ndolson/ont-guppy-cpu/bin
export PATH=$PATH:${BASEPATH}/msalit/ndolson/nanopore-pipeline

## Moving to scratch
cd ${BASEPATH}/msalit/ndolson

## activating env
source nanopore_env/bin/activate

## Creating run directory
JOBID=$(sacct -n -X -s r --format jobid --name ONT-test | sed 's/ //g')
mkdir -p ONT-pipe-run-logs/${JOBID}
cd ONT-pipe-run-logs/${JOBID}


## running running pipeline 
python ${BASEPATH}/msalit/ndolson/nanopore-pipeline/main.py $1

## deactivating env
deactivate
