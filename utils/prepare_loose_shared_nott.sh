#!/usr/bin
## usage prepare_loose_shared.sh INPUT.tar RUNNAME
## Prepared fast5s for pipeline from shared tars
#### Initial developed for Notts_GM24385 runs

echo "Creating fast5 tar balls"
tar_dir=$1
run_name=$2


## Compress Read directories
work_dir=/scratch/groups/msalit/nanopore/raw
cd ${tar_dir}/

run_dir=$(pwd)


cd ${run_dir}/reads/

## Looping through subread directories e.g. 0,1,2,...
for dir in */;
do
    cd ${run_dir}/reads/

    dir=${dir%*/}      # remove the trailing "/"
    echo ${read_name} ${dir}

    ## Grouping fast5s into 1000 read tar balls
        ls ${dir} > fast5_${dir}.txt
    lst_dir=$(pwd)
        split -l 1000 -d --additional-suffix=_fast5.lst \
            fast5_${dir}.txt ${run_name}_${dir}_
    rm fast5_${dir}.txt        
    echo ${lst_dir}
    ls -lh

    for lst in *_fast5.lst;
        do
            echo ${run_name} ${lst}
            out_tar=${lst%_fast5.lst}.tar
            cd ${lst_dir}/${dir}
        pwd
            tar -c --ignore-failed-read \
            --file=${work_dir}/${run_name}/fast5/${out_tar} \
            -T ${lst_dir}/${lst}
        rm ${lst_dir}/${lst}
    done
        echo "Done first loop"
pwd
done

