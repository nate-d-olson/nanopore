#!/usr/bin/env bash

set -o pipefail

combine_fast5s () {
    args=""
    subdirs=`find ${input_dir} -maxdepth 1 -type d`
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

cat_summary_seq(){
    echo "Concatenating seq summary files..."
    fastq_dirs=`find /oak/stanford/groups/msalit/nspies/nanopore/raw -maxdepth 2 -name "fastq"`

    if [ -e ${sequence_summary} ]; then
        rm ${sequence_summary}
    fi

    count=0
    for fastq_dir in ${fastq_dirs}; do
        ((count++))
        echo "${count} - ${fastq_dir}"
        cat ${fastq_dir}/*.sequence_summary.txt >> ${sequence_summary}
    done   
}

cat_bams(){
    ## Combine BAMs hs37d5
    echo "Concatenating combined bams for hs37d5..."
    ls ${input_dir}/*/aln_hs37d5/*combined.sorted.bam | \
        samtools cat | \
        samtools sort -@4 -m 2G -O bam \
        --reference ${hs37d5_ref} \
        -o ${bam37}

    echo "Performing hs37d5 QC..."
    python3 quick_qc.py ${bam37} ${output}.hs37d5.qc.pdf | tee ${output}.hs37d5.qc.txt

    ## Combine BAMs GRCh38
    echo "Concatenating combined bams for GRCh38..."
    ls ${input_dir}/*/aln_GRCh38/*combined.sorted.bam | \
        samtools cat | \
        samtools sort -@4 -m 2G -O bam \
        --reference ${grch38_ref} \
        -o ${bam38}
    
    samtools index ${bam38}

    echo "Performing GRCh38 QC..."
    python3 quick_qc.py ${bam38} ${output}.GRCh38.qc.pdf | tee ${output}.GRCh38.qc.txt
}

validate_combined_bams(){

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

input_dir=/scratch/groups/msalit/nanopore/nanopore-test/raw
output_dir=/scratch/groups/msalit/nanopore/nanopore-test/release/${output}
mkdir -p ${output_dir}
output=$output_dir/$output

## FAST5s
combine_fast5s

## FASTQs
fastq=${output}.fastq.gz
sequence_summary=${output}.sequence_summary.txt

cat_fastqs
validate_combined_fastq

## BAM files
hs37d5_ref=/oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa
grch38_ref=/oak/stanford/groups/msalit/shared/genomes/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
bam37=${output}_hs37d5.bam
bam38=${output}_GRCh38.bam

cat_bams
validate_combined_bams

