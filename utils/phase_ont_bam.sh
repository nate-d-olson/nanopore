#!/bin/sh
#SBATCH -p msalit
#SBATCH -o whatshap_ONTHG2tag_GRCh38_%A_%a.out
#SBATCH -e whatshap_ONTHG2tag_GRCh38_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH -t 02-00:00:00
#SBATCH -A msalit
module load python/3.6.1 biology samtools

whatshapdev=/oak/stanford/groups/msalit/ndolson/whatshap/whatshap_venv/bin/whatshap
releasedir=/scratch/groups/msalit/nanopore/release/ultra-long-ont
releasename=ultra-long-ont

## Phasing hs37d5 using whatshap
# refid=hs37d5
# phased_bam=${releasedir}/${releasename}_${refid}_phased.bam
# ${whatshapdev} haplotag -o ${phased_bam} \
#     -r /oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa \
#     /oak/stanford/groups/msalit/shared/giab/hg19/HG002/RTG.hg19.10x.trio-whatshap.vcf.gz \
#     ${releasedir}/${releasename}_${refid}_reheader.bam
# samtools index ${phased_bam}


# Phasing GRCh38 using whatshap
refid=GRCh38
phased_bam=${releasedir}/${releasename}_${refid}_phased.bam
${whatshapdev} haplotag -o ${phased_bam} \
    -r /oak/stanford/groups/msalit/shared/genomes/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
    /oak/stanford/groups/msalit/shared/giab/hg38/HG002/NA24385.GRCh38.phased_variants.reheadered.vcf.gz \
    ${releasedir}/${releasename}_${refid}_reheader.bam
samtools index ${phased_bam}
