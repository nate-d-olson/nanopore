#!/usr/bin/bash
#SBATCH -n 8 # number of cores
#SBATCH --mem 60g # memory pool for all cores
#SBATCH --job-name=Ashley37
#SBATCH --time=8:00:00
#SBATCH --partition=msalit,normal
#SBATCH --mail-type=ALL

BASEPATH=/oak/stanford/groups

## Loading required modules
module load python/3.6.1 R biology samtools

export PATH=$PATH:${BASEPATH}/msalit/ndolson/minimap2/

## Mapping to hs37d5
minimap2 -t 8 -aL -z 600,200 -x map-ont \
  -R '@RG\tID:GM24385_Ashley_Run1\tPU:PAD11290\tPL:nanopore\tPM:PromethION\tLB:GM24385_Ashley_Run1\tDT:2018-10-25 00:00:00\tPG:guppy-v2.3.1-dna_r9.4.1_450bps_flipflop.cfg\tDS:Flowcell=FLO-PRO002,kit=PSK-LSK109\tSM:HG002' \
   /oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa \
   /scratch/groups/msalit/nanopore/raw/GM24385_Ashley_Run1/GM24385_Ashley_Run1.fastq.gz | \
   samtools sort -m 2G -@8 -O bam \
   --reference /oak/stanford/groups/msalit/shared/genomes/hg19/hs37d5.fa \
       > /scratch/groups/msalit/nanopore/raw/GM24385_Ashley_Run1/aln_hs37d5/GM24385_Ashley_Run1.combined.sorted.bam 

samtools index /scratch/groups/msalit/nanopore/raw/GM24385_Ashley_Run1/GM24385_Ashley_Run1_hs37d5.combined.sorted.bam 

## Mapping to GRCh38
# minimap2 -t 12 -aL -z 600,200 -x map-ont \
#     -R '@RG\tID:GM24385_Ashley_Run1\tPU:PAD11290\tPL:nanopore\tPM:PromethION\tLB:GM24385_Ashley_Run1\tDT:2018-10-25 00:00:00\tPG:guppy-v2.3.1-dna_r9.4.1_450bps_flipflop.cfg\tDS:Flowcell=FLO-PRO002,kit=PSK-LSK109\tSM:HG002' \
#     /oak/stanford/groups/msalit/shared/genomes/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
#     /scratch/groups/msalit/nanopore/raw/GM24385_Ashley_Run1/GM24385_Ashley_Run1.fastq.gz | \
#     samtools sort -m 2G -@12 -O bam \
#     --reference /oak/stanford/groups/msalit/shared/genomes/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
#         > GM24385_Ashley_Run1_GRCh38.combined.sorted.bam 

# samtools index /scratch/groups/msalit/nanopore/raw/GM24385_Ashley_Run1/GM24385_Ashley_Run1_GRCh38.combined.sorted.bam 
