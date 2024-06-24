source("00_setup.R")

## DEREPLICATE IDENTICAL READS ----

# read in objects from last script
filt_fwd <- readRDS(paste(rdataPath, "filt_fwd.RData", sep = "/"))
filt_rev <- readRDS(paste(rdataPath, "filt_rev.RData", sep = "/"))

sink(file = paste(summaryPath, "/dereplicate_summary.txt", sep = ""))

derep_fwd <- derepFastq(filt_fwd, verbose = TRUE)
derep_rev <- derepFastq(filt_rev, verbose = TRUE)
# for 'Not all provided files exist' error:
# make sure you remove files from filtered that have zero post-filt reads
# table(file.exists(filt_fwd))
# Name the derep-class objects by the sample names
names(derep_fwd) <- sample_id
names(derep_rev) <- sample_id

sink()

saveRDS(derep_fwd,
        file = paste(rdataPath, "derep_fwd.RData", sep = "/"))
saveRDS(derep_rev,
        file = paste(rdataPath, "derep_rev.RData", sep = "/"))

