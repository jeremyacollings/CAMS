### Analysis of Microbial Community Data

set.seed(6)

library(rstan)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(patchwork)
library(vegan)
library(reshape2)

### Generating Data -----

##### Bringing data in -----

asv <- read.delim("Data/illumina_2023_asv-table.txt")
colnames(asv)[1]='ASV_Seq'

# other data attached to soil samples
# also includes extract number 
meta <- read.csv("Data/sample_order.csv")
meta$sp <- sub("_.*", "", meta$id)
meta$type <- ifelse(grepl("mono", meta$id), 'mono', "rhizo")

CAMS_ASV <- asv[, c(1, which(grepl("filt_J", colnames(asv))))]
names(CAMS_ASV) <- c("ASV_Seq", gsub("filt_", "", colnames(CAMS_ASV[,2:length(CAMS_ASV)])))
names(CAMS_ASV) <- c("ASV_Seq", sub("_.*", "",  colnames(CAMS_ASV[,2:length(CAMS_ASV)])))

# funguild IDs
guilds <- read.delim("Data/ASvforFG.taxa.guilds.txt")
colnames(guilds)[1]='ASV_Seq'

CAMS_ASV2 <- merge(guilds, CAMS_ASV, by="ASV_Seq") 

#filter out any identities less than "Probable" or "Highly Probable"
CAMS_ASV2 <- CAMS_ASV2 %>%
  filter(confidenceRanking != "Possible")

##### Making phyloseq object -----

# OTU table
otus <- as.matrix(CAMS_ASV2[, 18:length(CAMS_ASV2)])
rownames(otus) = CAMS_ASV2$ASV_Seq
# Taxa Table
tax <- as.matrix(CAMS_ASV2[,c(2:8, 11:12, 15)])
rownames(tax) = CAMS_ASV2$ASV_Seq
# Sample Data
sams <- meta
rownames(sams) <- sams$extract_num
CAMS_ps <- phyloseq(otu_table(otus, taxa_are_rows = TRUE), tax_table(tax), sample_data(sams))

# Pruning the ASVs
min_seq <- 5
min_samp <- .05
CAMS_ps2 <- prune_taxa(genefilter_sample(CAMS_ps, filterfun_sample(function(x) x > min_seq), A=min_samp*phyloseq::nsamples(CAMS_ps)), CAMS_ps)
CAMS_ps3 <- transform_sample_counts(CAMS_ps2, function(x) 1E6 * x/sum(x))

##### Subsetting phyloseq object by trophic mode -----

patho_ps <- subset_taxa(CAMS_ps3, trophicMode %in% c("Pathotroph", 
                                                      "Pathotroph-Saprotroph", 
                                                      "Pathotroph-Symbiotroph", 
                                                      "Pathotroph-Saprotroph-Symbiotroph"))
sapro_ps <- subset_taxa(CAMS_ps3, trophicMode %in% c("Saprotroph", 
                                                      "Pathotroph-Saprotroph", 
                                                      "Saprotroph-Symbiotroph", 
                                                      "Pathotroph-Saprotroph-Symbiotroph"))
symbio_ps <- subset_taxa(CAMS_ps3, trophicMode %in% c("Symbiotroph", 
                                                       "Pathotroph-Symbiotroph", 
                                                       "Saprotroph-Symbiotroph", 
                                                       "Pathotroph-Saprotroph-Symbiotroph"))
### Ordinations and PERMANOVAs -----

##### Running NMDSs and Extracting Scores ----- 

total_scores <- get_scores(CAMS_ps3, meta)
patho_scores <- get_scores(patho_ps, meta)
sapro_scores <- get_scores(sapro_ps, meta)
symbio_scores <- get_scores(symbio_ps, meta)

##### PERMANOVAS -----

total_test <- adonis2(distance(CAMS_ps3, method = "bray") ~ sp + type + sp:type,
                      data = as(sample_data(CAMS_ps3), "data.frame"))
total_test

patho_test <- adonis2(distance(patho_ps, method = "bray") ~ sp + type + sp:type,
                      data = as(sample_data(patho_ps), "data.frame"))
patho_test
sapro_test <- adonis2(distance(sapro_ps, method = "bray") ~ sp + type + sp:type,
                      data = as(sample_data(sapro_ps), "data.frame"))
sapro_test
symbio_test <- adonis2(distance(symbio_ps, method = "bray") ~ sp + type + sp:type,
                       data = as(sample_data(symbio_ps), "data.frame"))
symbio_test

##### Ordinations -----

# Species + Type Plots
total_sp_type_centroids <- get_sp_type_centroids(total_scores)
patho_sp_type_centroids <- get_sp_type_centroids(patho_scores)
sapro_sp_type_centroids <- get_sp_type_centroids(sapro_scores)
symbio_sp_type_centroids <- get_sp_type_centroids(symbio_scores)

sp_type_tot_ord <- plot_sp_type(total_scores, total_sp_type_centroids)

sp_type_patho_ord <- plot_sp_type(patho_scores, patho_sp_type_centroids)

sp_type_sapro_ord <- plot_sp_type(sapro_scores, sapro_sp_type_centroids)

plot_sp_type(symbio_scores, symbio_sp_type_centroids) + 
  theme(legend.position = "none")
symbio_sp_type_centroids <- get_sp_type_centroids(symbio_scores)
sp_type_symbio_ord <- plot_sp_type(symbio_scores, symbio_sp_type_centroids)

# Species Plots
total_sp_centroids <- get_sp_centroids(total_scores)
patho_sp_centroids <- get_sp_centroids(patho_scores)
sapro_sp_centroids <- get_sp_centroids(sapro_scores)
symbio_sp_centroids <- get_sp_centroids(symbio_scores)


plot_sp(total_scores, total_sp_centroids)
plot_sp(patho_scores, patho_sp_centroids)
plot_sp(sapro_scores, sapro_sp_centroids)
plot_sp(symbio_scores, symbio_sp_centroids)

# Type Plots
total_type_centroids <- get_type_centroids(total_scores)
patho_type_centroids <- get_type_centroids(patho_scores)
sapro_type_centroids <- get_type_centroids(sapro_scores)
symbio_type_centroids <- get_type_centroids(symbio_scores)

plot_type(total_scores, total_type_centroids)
plot_type(patho_scores, patho_type_centroids)
plot_type(sapro_scores, sapro_type_centroids)
plot_type(symbio_scores, symbio_type_centroids)

### Distance Plots -----

# distance matrices
tot_dist <- get_dist(CAMS_ps3@otu_table, meta)
patho_dist <- get_dist(patho_ps@otu_table, meta)
sapro_dist <- get_dist(sapro_ps@otu_table, meta)
symbio_dist <- get_dist(symbio_ps@otu_table, meta)

# mean distances compiled into dataframes
dist.dat <- cbind.data.frame(rbind(get_mean_dist(tot_dist), 
                                   get_mean_dist(patho_dist), 
                                   get_mean_dist(sapro_dist),
                                   get_mean_dist(symbio_dist)), 
                             com = rep(c("tot", "patho", "sapro", "symbio"), each = nrow(get_mean_dist(patho_dist))))
dist.dat$ord <- c(order(dist.dat$mean[which(dist.dat$com=="tot")]), 
                  order(dist.dat$mean[which(dist.dat$com=="patho")]), 
                  order(dist.dat$mean[which(dist.dat$com=="sapro")]), 
                  order(dist.dat$mean[which(dist.dat$com=="symbio")]))

# plots of sample distances by trophic mode
dist_box <- ggplot(data = dist.dat, aes(x = com, y = mean, color = com)) + 
  geom_jitter() + geom_boxplot(fill = NA) + 
  theme_classic(base_size = 12) + 
  xlab("Trophic Mode") + ylab("Bray-Curtis Dissimilarity") + 
  scale_color_manual(values = c("#988F2A", "#B18FCF", "#E4572E", "#2C2C34")) + 
  scale_x_discrete(labels = c("Pathotrophs", "Saprotrophs", 
                              "Symbiotrophs", "Total")) + 
  theme(axis.text.x = element_text(angle = 350, vjust = -.15), 
        legend.position = "none")

# plots of mean distance for each species pair for each trophic mode
mean_dist <- ggplot(data = dist.dat, aes(x = pair, y = mean, color = com)) + 
  geom_point(size = 2) + ylab("Mean Distance") + xlab("Species Pair") + 
  scale_color_manual(name = "Trophic Community", 
                     values = c("#988F2A", "#B18FCF", "#E4572E", "#2C2C34"), 
                     labels = c("Pathotrophs", "Saprotroph", "Symbiotroph", "Total")) + 
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = "none")

# plots of ordinal distance for each species pair for each trophic mode
ordn_dist <- ggplot(data = dist.dat, aes(x = pair, y = ord, color = com)) + 
  geom_point(size = 2) + ylab("Ordinal Distance") + xlab("Species Pair") + 
  scale_color_manual(name = "Trophic Community", 
                     values = c("#988F2A", "#B18FCF", "#E4572E", "#2C2C34"), 
                     labels = c("Pathotrophs", "Saprotroph", "Symbiotroph", "Total")) + 
  theme_classic(base_size = 12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = "none")

##### Estimating mean distance per pair -----

# only run if you haven't already generated these RDSs
sp_list <- cbind.data.frame(num = 1:8, 
                            sp = c("AA", "CA", "CC", "CG", 
                                   "CP", "NI", "PC", "PN"))

pulls.tot <- get_dist_est(CAMS_ps3)
write_rds(pulls.tot, "RDS_files/estDist_tot.RDS")
pulls.patho <- get_dist_est(patho_ps)
write_rds(pulls.patho, "RDS_files/estDist_patho.RDS")
pulls.sapro <- get_dist_est(sapro_ps)
write_rds(pulls.sapro, "RDS_files/estDist_sapro.RDS")
pulls.symbio <- get_dist_est(symbio_ps)
write_rds(pulls.symbio, "RDS_files/estDist_symbio.RDS")

##### Mantel Tests ----- 
tot.dists <- array(data = as.numeric(unlist(t(pulls.tot[,1:64]))), dim = c(8,8,nrow(pulls.tot)))
patho.dists <- array(data = as.numeric(unlist(t(pulls.patho[,1:64]))), dim = c(8,8,nrow(pulls.patho)))
sapro.dists <- array(data = as.numeric(unlist(t(pulls.sapro[,1:64]))), dim = c(8,8,nrow(pulls.sapro)))
symbio.dists <- array(data = as.numeric(unlist(t(pulls.symbio[,1:64]))), dim = c(8,8,nrow(pulls.symbio)))

sig.tot.patho <- c()
sig.tot.sapro <- c()
sig.tot.symbio <- c()
sig.patho.sapro <- c()
sig.patho.symbio <- c()
sig.sapro.symbio <- c()
for(i in sample(nrow(pulls.tot), 500)){
  sig.tot.patho <- c(sig.tot.patho, as.numeric(mantel(tot.dists[,,i], patho.dists[,,i])[4]))
  sig.tot.sapro <- c(sig.tot.sapro, as.numeric(mantel(tot.dists[,,i], sapro.dists[,,i])[4]))
  sig.tot.symbio <- c(sig.tot.symbio, as.numeric(mantel(tot.dists[,,i], symbio.dists[,,i])[4]))
  sig.patho.sapro <- c(sig.patho.sapro, as.numeric(mantel(patho.dists[,,i], sapro.dists[,,i])[4]))
  sig.patho.symbio <- c(sig.patho.symbio, as.numeric(mantel(patho.dists[,,i], symbio.dists[,,i])[4]))
  sig.sapro.symbio <- c(sig.sapro.symbio, as.numeric(mantel(sapro.dists[,,i], symbio.dists[,,i])[4]))
}
sum(sig.tot.patho < 0.05)/length(sig.tot.patho)
sum(sig.tot.sapro < 0.05)/length(sig.tot.sapro)
sum(sig.tot.symbio < 0.05)/length(sig.tot.symbio)
sum(sig.patho.sapro < 0.05)/length(sig.patho.sapro)
sum(sig.patho.symbio < 0.05)/length(sig.patho.symbio)
sum(sig.sapro.symbio < 0.05)/length(sig.sapro.symbio)

