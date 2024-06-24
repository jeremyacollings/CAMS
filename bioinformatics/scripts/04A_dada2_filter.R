source("00_setup.R")

## FILTER OUT LOW QUALITY READS ----

# set filtering parameters
maxee_fwd = 2
maxee_rev = 2
filt_name = "maxee2"
filt_label = "maxEE(2,2)" # for output figure

# read in objects from last script
trim_fwd <- readRDS(paste(rdataPath, "trim_fwd.RData", sep = "/"))
trim_rev <- readRDS(paste(rdataPath, "trim_rev.RData", sep = "/"))

# create file path for filtered reads to be created in next step
filt_fwd <- file.path(filtPath, paste("filt", gsub("trim_","",basename(trim_fwd)), sep="_"))
filt_rev <- file.path(filtPath, paste("filt", gsub("trim_","",basename(trim_rev)), sep="_"))

# standard filtering
filtered_reads <- filterAndTrim(trim_fwd, filt_fwd, trim_rev, filt_rev, 
                                maxN = 0, maxEE = c(maxee_fwd, maxee_rev), truncQ = 2, 
                                minLen = 50, rm.phix = TRUE, compress = F, 
                                multithread = FALSE, verbose = TRUE)
write_filtered_file(filtered_reads, filt_name, filt_label)

# filter out files that were empty after filtering
filt_fwd <- file.path(filtPath, paste(list.files(filtPath, pattern = "_R1.fastq.gz")))
filt_rev <- file.path(filtPath, paste(list.files(filtPath, pattern = "_R2.fastq.gz")))

saveRDS(filt_fwd,
        file = paste(rdataPath, "filt_fwd.RData", sep = "/"))
saveRDS(filt_rev,
        file = paste(rdataPath, "filt_rev.RData", sep = "/"))




