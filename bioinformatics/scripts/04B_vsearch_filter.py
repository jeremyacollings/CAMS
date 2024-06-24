#!/usr/bin/env python
# coding: utf-8

from glob import glob
import os

# CHANGE TO YOUR RUN NAME
run_name = "illumina_2022_relaxed"

path = os.getcwd() + "/relaxed_primer"
path_trimmed = f'{path}/trimmed_{run_name}'

out_name = f"vsearch_{run_name}.sh"

make_vsearch = open(out_name, 'wt') # name of bash script to be created

for file in os.listdir(path_trimmed): # location and file name (w/ wildcards) of your forward R1 reads
    sample = file[5:]
    cmd = "vsearch --fastx_filter" + file + \
        "--fastq_minlen 50" + \
        "--fastq_max22 2" +\
        "--log vsearch_filter_summary.txt" +\
        "--fastaout vsearch_" + sample_id + ".fasta"
    make_vsearch.write(cmd+'\n'+'\n')

make_vsearch.close()

vsearch --fastx_filter INPUT --fastq_minlen 50 --log vsearch_filter_summary.txt --fastaout vsearch_test_out.fasta --fastq_maxee 2

vsearch --fastq_stats --log filename.stats.out
