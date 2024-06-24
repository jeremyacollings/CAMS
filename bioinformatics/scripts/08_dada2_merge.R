source("00_setup.R")

## MERGED FWD AND REV READS ----

# read in objects from last script
derep_fwd <- readRDS(paste(rdataPath, "derep_fwd.RData", sep = "/"))
derep_rev <- readRDS(paste(rdataPath, "derep_rev.RData", sep = "/"))

# read in objects from last script
dada_fwd <- readRDS(paste(rdataPath, "dada_fwd.RData", sep = "/"))
dada_rev <- readRDS(paste(rdataPath, "dada_rev.RData", sep = "/"))

# merge reads
sink(file = paste(summaryPath, "/merging_summary.txt", sep = ""))

mergers <- mergePairs(dada_fwd, derep_fwd, 
                      dada_rev, derep_rev, verbose=TRUE)

sink()

saveRDS(mergers,
        file = paste(rdataPath, "mergers.RData", sep = "/"))

