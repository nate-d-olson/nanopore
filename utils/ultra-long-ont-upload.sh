#!/usr/bin/bash
#SBATCH -n 2 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --job-name=asp-ncbi
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --partition=msalit
 
## Adding aspera path
export PATH=/home/users/ndolson/.aspera/cli/bin:$PATH
 
## Key file location
SRAKEY=/home/users/ndolson/sra-rsakey
 
## Remote location - where files are being sent to
REMOTE=asp-sra@upload.ncbi.nlm.nih.gov:GIAB-incoming/NIST
 
ascp -v --overwrite=diff -i ${SRAKEY} -QT -l1000m -L aspera_log ultra-long-ont ${REMOTE}
