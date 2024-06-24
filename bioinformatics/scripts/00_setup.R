## SETUP

# user-specific changes:
# _ location of cutadapt
# _ primer set
# _ name of the sequencing run

# set name of this bioinformatics run, must match run_name in HTStream python script
run_name = "illumina_2023"

# load/install libraries and dependencies
library(dada2)
library(Biostrings)
library(ShortRead)
library(tidyverse)

# set directory
# setwd(paste("/gpfs/projects/roylab/cdelevic/", run_name, "/scripts", sep = "")) # wd set in slurm

# set locations to existing data
# cutadapt = "/Users/carolyndelevich/opt/anaconda3/envs/cutadaptenv/bin/cutadapt" # where is cutadapt (shell) downloaded?
cutadapt = "/gpfs/home/cdelevic/.conda/envs/illumina_env/bin/cutadapt"
rawPath = "../raw" # where are the raw reads?
primerPath = "../primers" # where are the fasta files containing the primer sequences

# set future file paths
demuxPath = paste("../demux", run_name, sep = "_") # where are the demultiplexed (and post HTStream) files?
prefiltPath = paste("../prefiltered", run_name, sep = "_") # where will the pre-filtered files go? (ambig bases)
rawprefiltPath = paste("../rawprefiltered", run_name, sep = "_") # where will the raw pre-filtered files go? (ambig bases)
filtPath = paste("../filtered", run_name, sep = "_")

trimPath0 = paste("../trimmed0", run_name, sep = "_") # where will cutadapt trimmed primer files go?
dir.create(trimPath0) # first folder for first cutadapt pass
trimPath = paste("../trimmed", run_name, sep = "_") # second folder for second cutadapt pass
dir.create(trimPath)

summaryPath = paste("../summary", run_name, sep = "_")
if (!dir.exists(summaryPath)){
  dir.create(summaryPath)
}

rdataPath = paste("./rdata", run_name, sep = "_")
if (!dir.exists(rdataPath)){
  dir.create(rdataPath)
}

# define primer set
primer_fwd = "CTTGGTCATTTAGAGGAAGTAA"  # ITS1F
primer_rev = "GCTGCGTTCTTCATCGATGC" #ITS2

# locate taxonomy reference database
taxonomydbPath = paste("../taxonomy-db", sep = "_")
unite_ref <- paste(taxonomydbPath, "sh_general_release_dynamic_29.11.2022.fasta", sep = "/") # where is unite ref db located

# functions
# summarize % seq loss per sample, export summary as csv
write_filtered_file <- function(filterAndTrim_output, file_name, filter_label){
  out <- filterAndTrim_output %>% 
    as.data.frame() %>%
    mutate(reads.lost = reads.in - reads.out) %>%
    mutate(perc.reads.lost = round((reads.lost/reads.in)*100,2)) %>%
    rownames_to_column(var = 'file.name')
  
  write.csv(out, paste(summaryPath, "/", file_name,".csv", sep = ""), row.names=F)
  
  ggplot(data = out, aes(x = perc.reads.lost)) +
    geom_histogram() +
    xlab("Percent Reads Lost to Filtering") +
    annotate("text", x = 90, y = -1, label = filter_label) +
    coord_cartesian(xlim = c(0, 100), clip = "off")
  
  ggsave(last_plot(),
         filename = paste(summaryPath, "/", file_name, ".jpeg", sep = ""))
}
