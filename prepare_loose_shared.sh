#!/usr/bin
## usage prepare_loose_shared.sh INPUT.tar RUNNAME
## Prepared fast5s for pipeline from shared tars
#### Initial developed for GM24385_B2_run#s
GETINFO=/oak/stanford/groups/msalit/ndolson/nanopore-pipeline/utils/get_flowcell_info.py
echo "Decompressing tar"

input_tar=$1
# tar -xf ${input_tar} 

echo "Getting flowcell info"
run_name=$2
fast5_file=$(find ${run_name}/ -name *fast5 -print -quit)
python ${GETINFO} ${fast5_file}

## Compress Read directories
mkdir fast5
cd ${run_name}/*/reads/
for dir in */;
do
    dir=${dir%*/}      # remove the trailing "/"
    echo ${dir}
    ls ${dir} > fast5_list.txt
    split -l 1000 -d --additional-suffix=_fast5.lst fast5_list.txt ${run_name}_${dir}_
    for lst in *.lst;
    do
        echo ${lst}
        out_tar=${lst%_fast5.lst}.tar ## replace _fast5.txt with tar
        cd ${dir}
        tar -cf ../../../../fast5/${out_tar} -T ../${lst}
        cd ../
    done
done
