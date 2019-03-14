#!/usr/bin
## usage prepare_loose_shared.sh INPUT.tar RUNNAME
## Prepared fast5s for pipeline from shared tars
#### Initial developed for GM24385_B2_run#s
GETINFO=/oak/stanford/groups/msalit/ndolson/nanopore-pipeline/utils/get_flowcell_info.py
echo "Decompressing tar"

input_tar=$1
tar -xf ${input_tar} 

echo "Getting flowcell info"
run_name=$2
fast5_file=$(ls ${run_name}/*/reads/0/*fast5 | head -1)
python ${GETINFO} ${fast5_file}

## Compress Read directories
cd ${run_name}/reads/*/
for dir in ${run_name}/reads/*/
do
    dir=${dir%*/}      # remove the trailing "/"
    tar -cf -a ../../fast5/${run_name}_${dir}.tar
done
