#!/usr/bin/env python
# coding: utf-8

from glob import glob
import os

# CHANGE TO YOUR RUN NAME
run_name = "illumina_2022"

path = os.getcwd()
path_raw = f'../raw' # raw folder should not include run name
path_trimmed = f'../trimmed_{run_name}'
path_stats = f'../stats_{run_name}'
path_demux = f'../demux_{run_name}'
path_primers = f'../primers'

out_name = f"./02_htstream_run.sh"

os.chdir(path)

if not os.path.exists(path_demux):
    # Create a new directory for trimmed files if it does not exist 
    os.makedirs(path_demux)
    print("Demultiplexed folder created.")
else:
    print("Demultiplexed folder already exists.")
    
if not os.path.exists(path_stats):
    # Create a new directory for stats files if it does not exist 
    os.makedirs(path_stats)
    print("Stats folder created.")
else:
    print("Stats folder already exists.")

make_hts = open(out_name, 'wt') # name of bash script to be created

for r1 in glob((path_raw + "/*_R1_*.gz")): # location and file name (w/ wildcards) of your forward R1 reads
    r2 = r1.replace("R1", "R2") # replicates R2 rev read file names based on R1 file name
    s = r1.split('/')[-1].replace("_L001_R1_001.fastq.gz", '') # extract sample ID
    cmd = "hts_Primers -1 " + r1 + " -2 " + r2 + " -f " + path_demux + "/demux_" + s + \
          " -P " + path_primers + "/climush_its_fwd.fa -Q " + path_primers + "climush_its_rev.fa \
          -l 5 -x -e 1 -d 6 -k -L " + path_stats + "/" + s + "_stats.txt"
    make_hts.write(cmd+'\n'+'\n')

make_hts.close()

