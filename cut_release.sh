#!/usr/bin/env bash

set -o pipefail

combine_fast5s () {
    args=""
    subdirs=`find ${input_dir} -maxdepth 1 -type d`
    # subdirs=`find ${input_dir} -maxdepth 1 -type d -printf "%T@ %f\n" | sort | cut -d' ' -f2`
    total_size=0
    chunk=1
    for run_dir in ${subdirs}; do
        [ ! -d ${run_dir}/fast5/ ] && continue
        # check to see if the current dir would put us past ~50GB
        cur_size=`du -s ${run_dir}/fast5/`
        cur_size=($cur_size)

        # if so, let's flush the fast5s so far
        echo ">>>>>>> ${run_dir} ${total_size}"
        if (( total_size+cur_size > 50000000 )); then
            echo "FLUSHING..."
            echo ${args}

            combined_fast5=${output}.raw_fast5s.${chunk}.tar
            python3 fast5_archives.py ${combined_fast5} ${args}
            let chunk++
            args=""
            total_size=0
        fi

        # add the current fast5 archives
        let total_size=${total_size}+cur_size
        for fast5 in `find ${run_dir}/fast5/ -name "*.tar"`; do
            cur_name=`basename ${run_dir}`
            args="${fast5},${cur_name} ${args}"
        done
    done

    # do one final flush
    echo "FINAL FLUSH..."
    echo ${args}
    python3 fast5_archives.py ${combined_fast5} ${args}
}


cat_fastqs () {
    echo "Concatenating fastqs..."
    fastq_dirs=`find /oak/stanford/groups/msalit/nspies/nanopore/raw -maxdepth 2 -name "fastq"`

    if [ -e ${fastq} ]; then
        rm ${fastq}
    fi

    
    count=0
    for fastq_dir in ${fastq_dirs}; do
        ((count++))
        echo "${count} - ${fastq_dir}"
        cat ${fastq_dir}/*.fastq.gz >> ${fastq}

        # if (( count > 1 )); then
        #     break
        # fi
    done
}

validate_combined_fastq () {
    echo "Validating fastq"
    if ! count=$(zcat ${fastq} | wc -l); then
        echo "Error in concatenated file!"
        exit $?
    fi

    read_count=$(( count / 4 ))
    echo "Collected ${read_count} reads"
}

map_reads () {
    declare -A assemblies=( ["hs37d5"]="/oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa"
        ["hg38"]="/oak/stanford/groups/msalit/shared/genomes/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna" )

    for assembly in "${!assemblies[@]}"; do
        echo "Mapping against ${assembly} - ${assemblies[$assembly]}"
        genome_path=${assemblies[$assembly]}
        out_cram=${output}.${assembly}.sorted.cram

        minimap2 -t ${threads} -a -z 600,200 -x map-ont ${genome_path} ${fastq}  \
           | samtools sort -m 1G -@${threads} -O cram --reference ${genome_path} > ${out_cram}
        samtools index ${out_cram}

        echo "Performing QC..."
        python3 quick_qc.py ${out_cram} ${output}.${assembly}.qc.pdf | tee ${output}.${assembly}.qc.txt
    done
}

if ! python3 -c "import numpy,tqdm"; then
    echo "Error: Need to run python >=3 with numpy installed"
    exit
fi


output=$1
threads=12

if [ -z ${output} ]; then
    echo "Need to specify the release name as argument 1 to this script"
    exit
fi

if [ "${output:0:1}" = "/" ]; then
    echo "Argument 1 should be a release name, not a full file path"
    exit
fi

input_dir=/oak/stanford/groups/msalit/nspies/nanopore/raw
output_dir=/oak/stanford/groups/msalit/nspies/nanopore/release/${output}
mkdir -p ${output_dir}
output=$output_dir/$output
fastq=${output}.fastq.gz

combine_fast5s
# cat_fastqs
# validate_combined_fastq
# map_reads