#!/usr/bin/bash
#SBATCH -n 2 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --job-name=qc-notebook
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

## running run qc notebook script
bash run_qc_notebook.sh 2018-07-17_qc.html

## deactivating env
deactivate
