#!/usr/bin/bash
#SBATCH -n 12 # number of cores
#SBATCH --mem 48g # memory pool for all cores
#SBATCH --job-name=cut_release
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL

## Loading required modules
module load python/3.6.1 R biology samtools

## Moving to scratch
cd /scratch/PI/msalit/ndolson

## activating env
source nspies_nanopore/bin/activate

## moving to nanopore repo
cd nanopore

## running run release script
bash cut_release.sh 2018-07-06

## deactivating env
deactivate
