source("00_setup.R")

## GET SAMPLE NAMES ----

# read in objects from last script
filt_fwd <- readRDS(paste(rdataPath, "filt_fwd.RData", sep = "/"))
filt_rev <- readRDS(paste(rdataPath, "filt_rev.RData", sep = "/"))

# extract sample names from trimmed files, if they have the format
sample_str_fwd <- as.vector(list.files(filtPath, pattern = "_R1.fastq.gz"), mode='any')
sample_id <- c() # for ASV table?
file_names_fwd <- character() # to create fwd filtered reads file name
for (i in 1:length(sample_str_fwd)) {
  spl_id <- strsplit(sample_str_fwd[i],"_")[[1]][2:6]
  file_names_fwd[i] <- substring(sample_str_fwd[i],6) # remove trim prefix
  sample_id[i] <- paste(spl_id, collapse = "_")
}

# to create rev filtered reads file name
sample_str_rev <- as.vector(list.files(filtPath, pattern = "_R2.fastq.gz"),mode='any')
file_names_rev <- character() # to create filtered reads file name
for (i in 1:length(sample_str_rev)) {
  file_names_rev[i] <- substring(sample_str_rev[i],6) # remove trim prefix
}

## ANALYZE QUALITY READ PROFILES ----

# plot any set of forward read profiles (n = 12)
plotQualityProfile(filt_fwd[1:12])
ggsave(paste(summaryPath, "qual-read-prof_fwd_sample.png", sep = "/"),
       last_plot(),
       width = 12,
       height = 8,
       units = "in")
# aggregate quality profile for all forward reads
plotQualityProfile(filt_fwd, aggregate = TRUE)
ggsave(paste(summaryPath, "qual-read-prof_fwd_aggr.png", sep = "/"),
       last_plot())

# plot any set of two reverse read profiles (n = 12)
plotQualityProfile(filt_rev[1:12])
ggsave(paste(summaryPath, "qual-read-prof_rev_sample.png", sep = "/"),
       last_plot(),
       width = 12,
       height = 8,
       units = "in")
# aggregate quality profile for all reverse reads
plotQualityProfile(filt_rev, aggregate = TRUE)
ggsave(paste(summaryPath, "qual-read-prof_rev_aggr.png", sep = "/"),
       last_plot())

## ERROR RATES ----

set.seed(111) # error rates learned from subset of data

error_fwd <- learnErrors(filt_fwd, nbases=1e8, multithread=TRUE)
# 103657569 total bases in 450953 reads 
# from 12 samples will be used for learning the error rates.
error_rev <- learnErrors(filt_rev, nbases=1e8, multithread=TRUE)
# 103734457 total bases in 450953 reads 
# from 12 samples will be used for learning the error rates.

# plot error rates to check
# red line - if quality scores exactly matched their nominal definition
plotErrors(error_fwd, nominalQ = TRUE)
ggsave(paste(summaryPath, "error-rates_fwd.png", sep = "/"),
       plot = last_plot())
plotErrors(error_rev, nominalQ = TRUE)
ggsave(paste(summaryPath, "error-rates_rev.png", sep = "/"),
       plot = last_plot())

saveRDS(error_fwd,
     file = paste(rdataPath, "error_fwd.RData", sep = "/"))
saveRDS(error_rev,
     file = paste(rdataPath, "error_rev.RData", sep = "/"))
