#!/bin/sh
#SBATCH -p msalit
#SBATCH -o whatshap_ONTHG2tag_%A_%a.out
#SBATCH -e whatshap_ONTHG2tag_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH -t 00-50:00:00
#SBATCH -A msalit
module load python/3.6.1 biology samtools

whatshapdev=/oak/stanford/groups/msalit/ndolson/whatshap/whatshap_venv/bin/whatshap
releasedir=/scratch/groups/msalit/nanopore/release/ultra-long-ont
releasename=ultra-long-ont
refid=hs37d5
phased_bam=${releasedir}/${releasename}_${refid}_phased.bam
#phased_bam=${releasedir}/test_${refid}_phased.bam

## Phasing using whatshap

${whatshapdev} haplotag -o ${phased_bam} \
    -r /oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa \
    /oak/stanford/groups/msalit/shared/giab/hg19/HG002/RTG.hg19.10x.trio-whatshap.vcf.gz \
    ${releasedir}/${releasename}_${refid}.bam

samtools index ${phased_bam}
