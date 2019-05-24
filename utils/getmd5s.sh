#!/usr/bin/sh
#SBATCH -n 1 # number of cores
#SBATCH --mem 4g # memory pool for all cores
#SBATCH --job-name=ont-md5s
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --partition=msalit

md5sum ultra-long-ont/* > md5sum.chk 
