#!/usr/bin/bash
#SBATCH -n 4 # number of cores
#SBATCH --mem 18g # memory pool for all cores
#SBATCH --job-name=rel_header
#SBATCH --time=1-00:00:00
#SBATCH --partition=msalit
#SBATCH --mail-type=ALL


## Loading required modules
module load biology samtools

cd /scratch/groups/msalit/nanopore/release/ultra-long-ont

fix_header(){
    BAM=$1
    HEADERFIX=${BAM%.bam}_header.sam
    OUTBAM=${BAM%.bam}_reheader.bam

    ## Replacing header
    samtools view -H ${BAM} | sed "s/guppy-v2.3.1/guppy-v2.3.5/" > ${HEADERFIX}

    samtools reheader -P ${HEADERFIX} ${BAM} > ${OUTBAM}

    ## Indexing bam with fixed header
    samtools index ${OUTBAM}

    ## Sanity check
    samtools quickcheck ${OUTBAM}
}

## Fixing header for 37
fix_header ultra-long-ont_hs37d5.bam

## Fixing header for 38
fix_header ultra-long-ont_GRCh38.bam
