#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --job-name=QC-stats
#SBATCH --time=1-00:00:00
#SBATCH --partition=msalit,normal
#SBATCH --mail-type=ALL

module load python/3.6.1 R biology samtools
BASEPATH=/oak/stanford/groups
export PATH=$PATH:${BASEPATH}/msalit/ndolson/minimap2/
export PATH=$PATH:${BASEPATH}/msalit/ndolson/pigz/
export PATH=$PATH:${BASEPATH}/msalit/ndolson/ont-guppy-cpu/bin
export PATH=$PATH:${BASEPATH}/msalit/ndolson/nanopore-pipeline/

## Moving to scratch
cd ${BASEPATH}/msalit/ndolson

## activating env
source nanopore_env/bin/activate

cd /scratch/groups/msalit/nanopore/release/ultra-long-ont

python3 /oak/stanford/groups/msalit/ndolson/nanopore-pipeline/quick_qc.py \
		ultra-long-ont_GRCh38.bam \
		ultra-long-ont_GRCh38.qc.pdf | tee ultra-long-ont.GRCh38.qc.txt

python3 /oak/stanford/groups/msalit/ndolson/nanopore-pipeline/quick_qc.py \
		ultra-long-ont_GRCh38.bam \
		ultra-long-ont_GRCh38.qc.pdf | tee ultra-long-ont.GRCh38.qc.txt

## deactivating env
deactivate