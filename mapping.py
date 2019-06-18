import subprocess
import sys


def run_command(command):
    sys.stderr.write(f"Running command '{command}'...")
    sys.stderr.flush()
    
    subprocess.check_call(command, shell=True)


def run_mapping(fastq, out_bam, read_group_args, genome_path, threads=1):
    minimap2 = "minimap2"

    read_group = (
        f"@RG\\tID:{read_group_args['run']}\\t"
        f"PU:{read_group_args['flowcell_id']}\\t"
        f"PL:nanopore\\t"
        f"PM:{read_group_args['platform_model']}\\t"
        f"LB:{read_group_args['sample']}\\t"
        f"DT:{read_group_args['date']}\\t"
        f"PG:guppy-v2.3.5-{read_group_args['guppy_config']}\\t"
        f"DS:Flowcell={read_group_args['flowcell_type']},kit={read_group_args['flowcell_kit']}\\t"
        f"SM:HG002"
    )

    # we specify -z because the default seems to spit out alignments that
    # are not particularly contiguous 
    map_command = f"{minimap2} -t {threads} -aL -z 600,200 -x map-ont -R \'{read_group}\' {genome_path} {fastq} " \
                  f"| samtools sort -m 1G -@{threads} -O bam --reference {genome_path} > {out_bam}\n"

    run_command(map_command)


    index_command = f"samtools index {out_bam}"
    run_command(index_command)

    ## BAM file sanity checks ####################
    ## Checking file - mapping stats and for incomplete file
    idxstats_command = f"samtools idxstats {out_bam}"
    run_command(idxstats_command)

    ## Checking header format
    header_command = f"samtools view -H {out_bam}"
    run_command(header_command) 


def merge_bams(combined_path, bams, genome_path):
    # merge_command = f"samtools merge -f -O cram --reference {genome_path} {combined_path} {' '.join(bams)}"
    # run_command(merge_command)

    merge_command = f"samtools cat {' '.join(bams)} | samtools sort -@4 -m 2G -O bam --reference {genome_path} -o {combined_path}"
    run_command(merge_command)
    
    run_command(f"samtools index {combined_path}")
