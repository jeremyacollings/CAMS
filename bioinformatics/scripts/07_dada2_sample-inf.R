source("00_setup.R")

## SAMPLE INFERENCE ----

# read in objects from last script
derep_fwd <- readRDS(paste(rdataPath, "derep_fwd.RData", sep = "/"))
derep_rev <- readRDS(paste(rdataPath, "derep_rev.RData", sep = "/"))

error_fwd <- readRDS(paste(rdataPath, "error_fwd.RData", sep = "/"))
error_rev <- readRDS(paste(rdataPath, "error_rev.RData", sep = "/"))

# dereplicate
dada_fwd <- dada(derep_fwd, err = error_fwd, multithread = TRUE)
dada_rev <- dada(derep_rev, err = error_rev, multithread = TRUE)

saveRDS(dada_fwd,
        file = paste(rdataPath, "dada_fwd.RData", sep = "/"))
saveRDS(dada_rev,
        file = paste(rdataPath, "dada_rev.RData", sep = "/"))