#!/bin/sh
#SBATCH -n 8 # number of cores
#SBATCH --mem 16g # memory pool for all cores
#SBATCH --job-name=guppy_sh
#SBATCH --time=10:00:00
#SBATCH --partition=msalit,owners
#SBATCH --mail-type=ALL

fast5path=/oak/stanford/groups/msalit/ndolson/nanopore-test/raw/gm24385_181205DNA2e/fast5  
outdir=/oak/stanford/groups/msalit/ndolson/nanopore-test/raw/gm24385_181205DNA2e/fastq
flowcell=FLO-MIN106
kit=SQK-RAD004
threads=12

## Run guppy on multi-seq fast5s 
singularity exec \
    /home/groups/msalit/ndolson/simg/guppy_latest.sif \
    guppy_basecaller \
    -c dna_r9.4.1_450bps_flipflop.cfg \
    --input_path ${fast5path} \
    --save_path ${outdir} \
    --num_callers ${threads}
    # --flowcell ${flowcell} \
    # --kit ${kit} \
