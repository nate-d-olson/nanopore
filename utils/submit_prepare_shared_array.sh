#!/usr/bin/bash
#SBATCH --job-name=ont-prep-array
#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
#SBATCH --array=1-9
#SBATCH --time=08:00:00
#SBATCH --partition=msalit,normal
#SBATCH -n 4 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL

# Print this sub-job's task ID
echo "Prepare shared array job: " $SLURM_ARRAY_TASK_ID

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

cd /oak/stanford/groups/msalit/nspies/nanopore/nanopore_shared/loman_lab

bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh \
	GIAB_GM24385_UB_Run${SLURM_ARRAY_TASK_ID}.tar \
	GIAB_GM24385_UB_Run${SLURM_ARRAY_TASK_ID}

deactivate
