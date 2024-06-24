source("00_setup.R")

# do you want to see the number of primer occurrences in raw reads? (yes if TRUE)
raw_read_primers = FALSE

# do you want to cut primers for paired-end or fwd/rev separate? (paired if TRUE)
cutadapt_pairend = FALSE

# extract file names of fwd and rev demux files in HTStream output folder
demux_fwd <- sort(list.files(demuxPath, pattern = "_R1.fastq.gz", full.names = TRUE))
demux_rev <- sort(list.files(demuxPath, pattern = "_R2.fastq.gz", full.names = TRUE))

# randomly decides which sample to show printout for primer table
# sample_num = sample.int(length(demux_fwd), size = 1) # doesn't seem to want to work
sample_num = 1

## PRE-FILTER ----

print("Running prefilter step...")

# if the prefilter folder is not empty, empty it before next step
if (length(list.files(path=prefiltPath)) != 0) {
  print("Prefilt folder not empty, emptying now...")
  for (file in list.files(path = prefiltPath)) {
    file.remove(paste(prefiltPath, file, sep = "/"))
  }
}

# get file path of each sample
# I updated this so that you can have any length run name
prefilt_fwd <- file.path(prefiltPath, paste("prefilt", gsub("demux_","",basename(demux_fwd)), sep = "_"))
prefilt_rev <- file.path(prefiltPath, paste("prefilt", gsub("demux_","",basename(demux_rev)), sep = "_"))
  
# pre-filter sequences for ambiguous bases (N)
filterAndTrim(demux_fwd, prefilt_fwd,
              demux_rev, prefilt_rev,
              maxN = 0, multithread = FALSE)

# CHECK FOR PRESENCE OF PRIMERS ----

print("Checking for presence of primers...")

# # get all possible orientations of fwd and rev primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(primer_fwd)
REV.orients <- allOrients(primer_rev)

# count number of times and orientation of primers appearing in seqs
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

sink(file = paste(summaryPath, "/primers_prefilt.txt", sep = ""))

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = prefilt_fwd[[sample_num]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = prefilt_rev[[sample_num]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = prefilt_fwd[[sample_num]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = prefilt_rev[[sample_num]]))

sink()

## COMPARE TO RAW READS (OPTIONAL) ----

if (raw_read_primers == TRUE){
  # if there are issues removing primers from seq reads, it can be helpful to
  # compare to primer presence/orientation in raw reads (pre-HTStream)
  raw_fwd <- sort(list.files(rawPath, pattern = "_R1_001.fastq.gz", full.names = TRUE))
  raw_rev <- sort(list.files(rawPath, pattern = "_R2_001.fastq.gz", full.names = TRUE))
  
  # pre-filter sequences for ambiguous bases (N)
  rawprefilt_fwd <- file.path(rawprefiltPath, paste("prefilt", basename(raw_fwd), sep = "_"))
  rawprefilt_rev <- file.path(rawprefiltPath, paste("prefilt", basename(raw_rev), sep = "_"))
  
  filterAndTrim(raw_fwd, rawprefilt_fwd,
                raw_rev, rawprefilt_rev,
                maxN = 0, multithread = FALSE)
  
  sink(file = paste(summaryPath, "/primers_raw.txt", sep = ""))
  
  rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = rawprefilt_rev[[sample_num]]),
        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = rawprefilt_rev[[sample_num]]),
        REV.ForwardReads = sapply(REV.orients, primerHits, fn = rawprefilt_rev[[sample_num]]),
        REV.ReverseReads = sapply(REV.orients, primerHits, fn = rawprefilt_rev[[sample_num]]))
  
  sink()
}

## REMOVE PRIMERS WITH CUTADAPT ----

# if HTStream doesn't properly remove primers, get rid of the primer removal
# step via HTStream (e.g., -k statement) and run cutadapt

print("Running cutadapt...")

# export the output of the cutadapt to a text file
sink(file = paste(summaryPath, "/cutadapt_summary.txt", sep = ""))

system2(cutadapt, args = '--version') # check cutadapt located

trim0_fwd <- file.path(trimPath0, paste("trim", gsub("demux_","",basename(demux_fwd)), sep = "_"))
trim0_rev <- file.path(trimPath0, paste("trim", gsub("demux_","",basename(demux_rev)), sep = "_"))

trim_fwd <- file.path(trimPath, paste("trim", gsub("demux_","",basename(demux_fwd)), sep = "_"))
trim_rev <- file.path(trimPath, paste("trim", gsub("demux_","",basename(demux_rev)), sep = "_"))

# get reverse complement of each fwd/rev primer
primer_fwd_rc <- dada2:::rc(primer_fwd)
primer_rev_rc <- dada2:::rc(primer_rev)

# cutadapt paired-end trimming
if (cutadapt_pairend == TRUE){
  fwd_totrim <- paste("-g", primer_fwd, "-a", primer_rev_rc)
  rev_totrim <- paste("-G", primer_rev, "-A", primer_fwd_rc)
  for (i in seq_along(demux_fwd)) {
    system2(cutadapt, args = c(fwd_totrim, rev_totrim, "-n", 2, # -n 2 required to remove FWD and REV from reads (runs check twice)
                               "-o", trim_fwd[i], "-p", trim_rev[i], # trimmed output files
                               prefilt_fwd[i], prefilt_rev[i])) # prefiltered input files
  }
}else{
  print("Starting cutadapt for forward and reverse reads separately...")
  
  # cutadapt forward read trimming - pass 01
  fwd_totrim <- paste("-g", primer_fwd, "-a", primer_rev_rc)
  for (i in seq_along(demux_fwd)) {
    system2(cutadapt, args = c(fwd_totrim, "-n", 2, # -n 2 required to remove FWD and REV from reads (runs check twice)
                               "-o", trim0_fwd[i], # trimmed output files
                               prefilt_fwd[i]))
  }
  
  # cutadapt reverse read trimming - pass 01
  rev_totrim <- paste("-g", primer_rev, "-a", primer_fwd_rc)
  for (i in seq_along(demux_rev)) {
    system2(cutadapt, args = c(rev_totrim, "-n", 2, # -n 2 required to remove FWD and REV from reads (runs check twice)
                               "-o", trim0_rev[i], # trimmed output files
                               prefilt_rev[i])) # prefiltered input files
  }
  
  # switch these out, then run cutadapt again
  temp1 <- primer_rev
  temp2 <- primer_fwd_rc
  primer_rev <- primer_fwd
  primer_fwd_rc <- primer_rev_rc
  primer_fwd <- temp1
  primer_rev_rc <- temp2
  
  # cutadapt forward read trimming - pass 02
  fwd_totrim <- paste("-g", primer_fwd, "-a", primer_rev_rc)
  for (i in seq_along(demux_fwd)) {
    system2(cutadapt, args = c(fwd_totrim, "-n", 2, # -n 2 required to remove FWD and REV from reads (runs check twice)
                               "-o", trim_fwd[i], # trimmed output files
                               trim0_fwd[i]))
  }
  
  # cutadapt reverse read trimming
  rev_totrim <- paste("-g", primer_rev, "-a", primer_fwd_rc)
  for (i in seq_along(demux_rev)) {
    system2(cutadapt, args = c(rev_totrim, "-n", 2, # -n 2 required to remove FWD and REV from reads (runs check twice)
                               "-o", trim_rev[i], # trimmed output files
                               trim0_rev[i])) # prefiltered input files
  }
}

# check whether primers were properly removed - pass 01 (should not be all zeros)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = trim0_fwd[[sample_num]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = trim0_rev[[sample_num]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = trim0_fwd[[sample_num]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = trim0_rev[[sample_num]]))

# check whether primers were properly removed - pass 02 (should be all zeros)
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = trim_fwd[[sample_num]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = trim_rev[[sample_num]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = trim_fwd[[sample_num]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = trim_rev[[sample_num]]))

sink()

# save objects required for the next script
saveRDS(trim_fwd,
        file = paste(rdataPath, "trim_fwd.RData", sep = "/"))
saveRDS(trim_rev,
        file = paste(rdataPath, "trim_rev.RData", sep = "/"))
