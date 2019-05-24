#!/usr/bin/sh
python3 /oak/stanford/groups/msalit/ndolson/nanopore-pipeline/quick_qc.py \
		ultra-long-ont_GRCh38.bam \
		ultra-long-ont_GRCh38.qc.pdf | tee ultra-long-ont.GRCh38.qc.txt
