#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --job-name=ONT
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL

BASEPATH=/oak/stanford/groups

## Loading required modules
module load python/3.6.1 R biology samtools

export PATH=$PATH:${BASEPATH}/msalit/ndolson/minimap2/:${BASEPATH}/msalit/ndolson/pigz/

## Moving to scratch
cd ${BASEPATH}/msalit/ndolson

## activating env
source nspies_nanopore/bin/activate

## moving to nanopore repo
cd nanopore

## running running pipeline 
python main.py $1

## deactivating env
deactivate
