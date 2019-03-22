#!/usr/bin
## usage prepare_loose_shared.sh INPUT.tar RUNNAME
## Prepared fast5s for pipeline from shared tars
#### Initial developed for GM24385_B2_run#s
GETINFO=/oak/stanford/groups/msalit/ndolson/nanopore-pipeline/utils/get_flowcell_info.py
echo "Decompressing tar"

input_tar=$1

## Tar -z for gzipped tar balls
if [ "$input_tar" == "*.tar.gz" ]
then
	tar -xzf ${input_tar} 
else
	tar -xf ${input_tar}
fi

echo "Getting flowcell info"
run_name=$2
fast5_file=$(find -name *fast5 ${run_name} -print -quit)

echo ${fast5_file}

python ${GETINFO} ${fast5_file}

## Compress Read directories
mkdir ${run_name}_fast5
work_dir=$(pwd)
cd ${run_name}/

## Looping through read directories GA000#
run_dir=$(pwd)

for read_dir in */;
do
    read_dir=${read_dir%*/}
    cd ${run_dir}/${read_dir}/reads/

    ## Looping through subread directories e.g. 0,1,2,...
    for dir in */;
    do
        cd ${run_dir}/${read_dir}/reads/

        dir=${dir%*/}      # remove the trailing "/"
        echo ${read_dir} ${dir}
	
	## Grouping fast5s into 1000 read tar balls
        ls ${dir} > fast5_${dir}.txt
	lst_dir=$(pwd)
        split -l 1000 -d --additional-suffix=_fast5.lst fast5_${dir}.txt ${run_name}_${read_dir}_${dir}_
	rm fast5_${dir}.txt        
	echo ${lst_dir}
	ls -lh
        
	for lst in *_fast5.lst;
    	do
            echo ${read_dir} ${lst}
            out_tar=${lst%_fast5.lst}.tar
            cd ${lst_dir}/${dir}
	    pwd
            tar -c --ignore-failed-read --file=${work_dir}/${run_name}_fast5/${out_tar} -T ${lst_dir}/${lst}
  	    rm ${lst_dir}/${lst}
	done
        echo "Done first loop"
	pwd
    done
    echo "Done second loop"
    pwd
done
