#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 12g # memory pool for all cores
#SBATCH --job-name=Prep-Loose
#SBATCH --time=1-00:00:00
#SBATCH --partition=msalit,normal
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

cd /oak/stanford/groups/msalit/nspies/nanopore/nanopore_shared/loman_lab

#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run1.tar GIAB_GM24385_UB_Run1
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run2.tar GIAB_GM24385_UB_Run2
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run3.tar GIAB_GM24385_UB_Run3
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run4.tar GIAB_GM24385_UB_Run4
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run5.tar GIAB_GM24385_UB_Run5
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run6.tar GIAB_GM24385_UB_Run6
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run7.tar GIAB_GM24385_UB_Run7
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run8.tar GIAB_GM24385_UB_Run8
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run9.tar GIAB_GM24385_UB_Run9
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run10.tar GIAB_GM24385_UB_Run10
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run11.tar GIAB_GM24385_UB_Run11
#bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GIAB_GM24385_UB_Run12.tar GIAB_GM24385_UB_Run12

## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run1.tar.gz GM24385_B2_run1
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run2.tar.gz GM24385_B2_run2
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run3.tar.gz GM24385_B2_run3
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run3_4.tar.gz GM24385_B2_run3_4
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run4.tar.gz GM24385_B2_run4
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run5.tar.gz GM24385_B2_run5
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run6.tar.gz GM24385_B2_run6
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run7.tar.gz GM24385_B2_run7
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run8.tar.gz GM24385_B2_run8
## In progress: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run9.tar.gz GM24385_B2_run9
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run10.tar.gz GM24385_B2_run10
## Sync In progress: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run11.tar.gz GM24385_B2_run11
## Completed: bash ${BASEPATH}/msalit/ndolson/nanopore-pipeline/utils/prepare_loose_shared.sh GM24385_B2_run12.tar.gz GM24385_B2_run12

deactivate
