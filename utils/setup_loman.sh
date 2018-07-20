#!/usr/bin/bash
BASEDIR=/oak/stanford/groups/msalit/nspies/nanopore
DATROOT=GIAB_GM24385_UB
for i in {1..6} 8;
  do
    RUNNAME=${DATROOT}_Run${i}
    echo "Setting up ${RUNNAME}"
	
    ln -sf ${BASEDIR}/nanopore_shared/loman_lab/${RUNNAME}.tar \
      ${RUNNAME}/fast5

    cd ${RUNNAME}/fast5

    ## Get list of fastq files
    echo "Extracting fastq" 
    tar -v --extract \
      --file=${RUNNAME}.tar \
      --wildcards --no-anchored '*fastq' '*summary*.txt'
    

    echo "Moving and compressing fastq files"
    for j in ${RUNNAME}/*/*{fastq,txt}; 
	do
		NEWNAME=$(echo ${j} | sed 's/${RUNNAME}\///')
		NEWNAME=$(echo ${NEWNAME} | sed 's/\//_/g')
    		mv ${j} ../fastq/${NEWNAME} 
	done
    rm -r ${RUNNAME}



    echo "Seq Info"
    cd ../fastq
    head *sequencing_summary_0.txt


    ## Compress fastq
    gzip *fastq

    ## Moving back to raw data directory
    cd ${BASEDIR}/raw
  done
