source("00_setup.R")

## ASSIGN TAXONOMY ----

# load in objects from previous scripts
seq_table_nochim <- readRDS(paste(rdataPath, "seq_table_nochim.RData", sep = "/"))

# the unite repository needs to be downloaded from unite and saved locally
# won't recognize _dev version
if(length(list.files(taxonomydbPath)) == 0) 
  stop("Error: Taxonomic database not downloaded to taxonomy-db folder.")

taxa <- assignTaxonomy(seq_table_nochim, unite_ref, 
                       multithread = TRUE, tryRC = TRUE) #tryRC=try rev complement

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tax_fname = paste(summaryPath, "/", run_name, "_taxonomy.txt", sep = "")
write.table(taxa, tax_fname, 
            sep='\t', row.names = TRUE, quote = FALSE) # write the taxa table out to a file just in case -can delete later

## CREATE ASV FILE ----

seq_table_nochim.t = t(seq_table_nochim)
asv_fname = paste(summaryPath, "/", run_name, "asvs.txt", sep = "")
write.table(seq_table_nochim.t, asv_fname,
            sep="\t", row.names=TRUE, col.names=NA, quote=FALSE) #write the ASV table out to a file just in case

## GET REPRESENTATIVE SEQUENCES ----

# we can also get the representative sequences in case we want to do phylogenetic analyses
uniquesToFasta(seq_table_nochim, 
               fout = paste(summaryPath, "/", run_name, "rep-seqs.fna", sep = ""), 
               ids=colnames(seq_table_nochim))

## CONSTRUCT ASV TABLE ----

asv <- read.delim(asv_fname)
taxa <- read.delim(tax_fname)

# reassign seqs to X in taxa table for merge
taxa$X <- rownames(taxa)

asv_table <- merge(taxa, asv, by = "X", sort = FALSE)

# do some basic filtering of the table and save
asvtab_fname = paste(summaryPath, "/", run_name, "_asv-table.txt", sep = "")
asv_table %>%
  subset(Kingdom != "k__unidentified") %>% # remove unidentified
  subset(Kingdom != "k__Rhizaria") %>% # remove rhizaria ****for ITS only
  write.table(asvtab_fname, sep='\t', 
              row.names = FALSE, quote = FALSE)

## GET SAMPLE SUMS FOR RAREFACTION ----

sum0 <- data.frame(X = colSums(asv_table[9:ncol(asv_table)])) # get sample sums for rarefaction
sum1 <- arrange(sum0, desc(X))
rare_fname = paste(summaryPath, "/", run_name, "_rarefaction.txt", sep = "")
write.table(sum1, file = rare_fname, sep="\t") #change file name --sample sum table
