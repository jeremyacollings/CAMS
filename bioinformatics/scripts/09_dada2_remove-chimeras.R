source("00_setup.R")

## CONSTRUCT SEQUENCE TABLE ----

# read in objects from last script
mergers <- readRDS(paste(rdataPath, "mergers.RData", sep = "/"))

dada_fwd <- readRDS(paste(rdataPath, "dada_fwd.RData", sep = "/"))
dada_rev <- readRDS(paste(rdataPath, "dada_rev.RData", sep = "/"))

# make seq table
seq_table <- makeSequenceTable(mergers)
dim(seq_table)

## REMOVE CHIMERAS ----
sink(file = paste(summaryPath, "/remove-chimeras_summary.txt", sep = ""))

seq_table_nochim <- removeBimeraDenovo(seq_table, method = "consensus", 
                                       multithread = multithread_bool, verbose = TRUE)

sink()

saveRDS(seq_table_nochim,
        file = paste(rdataPath, "seq_table_nochim.RData", sep = "/"))

# check distribution of sequence lengths
seq_lengths <- table(nchar(getSequences(seq_table_nochim)))
write.csv(seq_lengths,
          file = paste(summaryPath, "seq-lengths.csv", sep = "/"),
          row.names = FALSE)


## TRACK READS THROUGH PIPELINE ----

getN <- function(x) sum(getUniques(x))
track_reads <- cbind(filtered_reads,
                     sapply(dada_fwd, getN), sapply(dada_rev, getN), 
                     sapply(mergers, getN),
                     rowSums(seq_table_nochim))
colnames(track_reads) <- c("input", "filtered", "denoised_fwd", 
                           "denoised_rev", "merged", "nonchim")
# this table will include the samples that were lost after filtering, check 
# sample names also include these samples
track_read_rn <- rownames(track_reads)
sample_id_track <- list()
for (i in 1:length(track_read_rn)) {
  spl_id <- strsplit(track_read_rn[i],"_")[[1]][2:6]
  sample_id_track[i] <- paste(spl_id, collapse = "_")
}
rownames(track_reads) <- sample_id_track

track_reads %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_id_fwd') %>%
  write.csv(file = paste(summaryPath, "pipeline-read-tracking.csv", sep = "/"),
            row.names = FALSE)

