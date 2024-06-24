
### Functions

##### Microbial Analyses -----

# run NMDS and extract scores
get_scores <- function(dat, samps){
  # dat - phyloseq object
  # samps - metadata
  m1 <- ordinate(dat, "NMDS", "bray")
  # extracting points
  mds_dat <- as.data.frame(m1$points)
  mds_dat$extract_num <- row.names(mds_dat)
  mds_dat <- merge(mds_dat, samps, by = "extract_num")
  mds_dat$sp.type <- paste(mds_dat$sp, mds_dat$type, sep = ".")
  mds_dat
}

# Calculate centroids by Species and Sample Type
get_sp_type_centroids <- function(scores){
  # scores - output of get_scores()
  sp.type.mean1 <- aggregate(scores$MDS1, list(scores$sp.type), mean)
  sp.type.mean2 <- aggregate(scores$MDS2, list(scores$sp.type), mean)
  sp.type.se1 <- aggregate(scores$MDS1, list(scores$sp.type), function(x) sd(x)/sqrt(length(x)))
  sp.type.se2 <- aggregate(scores$MDS2, list(scores$sp.type), function(x) sd(x)/sqrt(length(x)))
  sp.type.dat <- cbind.data.frame(sp.type = sp.type.mean1[1], 
                                  sp = rep(sort(unique(scores$sp)), each = 2), 
                                  type = rep(c("mono", "rhizo"), 8), 
                                  mean1 = sp.type.mean1[2], 
                                  mean2 = sp.type.mean2[2], 
                                  se1 = sp.type.se1[2], 
                                  se2 = sp.type.se2[2])
  colnames(sp.type.dat) <- c("sp.type", "sp", "type", "mean1", "mean2", "se1", "se2")
  sp.type.dat
}

# Plot NMDS in ggplot by Species and Sample Type
plot_sp_type <- function(scores, centroids){
  # scores - output from get_scores()
  # centroids - output from get_sp_type_centroids()
  ggplot() + 
    geom_point(data = scores, aes(x = MDS1, y = MDS2, 
                                  color = sp, shape = type), alpha = .6) + 
    geom_point(data = centroids, aes(x = mean1, y = mean2,
                                     color = sp, shape = type), size = 3) + 
    geom_errorbar(data = centroids, aes(x = mean1, y = mean2, 
                                        xmin = mean1 - se1, xmax = mean1 + se1,
                                        color = sp), width = 0) + 
    geom_errorbar(data = centroids, aes(x = mean1, y = mean2, 
                                        ymin = mean2 - se2, ymax = mean2 + se2,
                                        color = sp), width = 0) +
    theme_classic(base_size = 14) + 
    scale_shape_manual(name = "Sample Type", values = c(16, 17), 
                       labels = c("Monoculture Plot", "Experimental Plot")) + 
    scale_color_manual(name = "Species", 
                       values = c("#4381A4", "#48A9A6", "#E4DFDA", "#D4B484", 
                                  "#C1666B", "#896F69","#F09021", "#AA1E57"))
}

# Calculate centroids by Species
get_sp_centroids <- function(scores){
  # scores - output from get_scores
  sp.mean1 <- aggregate(scores$MDS1, list(scores$sp), mean)
  sp.mean2 <- aggregate(scores$MDS2, list(scores$sp), mean)
  sp.se1 <- aggregate(scores$MDS1, list(scores$sp), function(x) sd(x)/sqrt(length(x)))
  sp.se2 <- aggregate(scores$MDS2, list(scores$sp), function(x) sd(x)/sqrt(length(x)))
  sp.dat <- cbind.data.frame(sp = sp.mean1[1],  
                             mean1 = sp.mean1[2], 
                             mean2 = sp.mean2[2], 
                             se1 = sp.se1[2], 
                             se2 = sp.se2[2])
  colnames(sp.dat) <- c("sp", "mean1", "mean2", "se1", "se2")
  sp.dat
}

# Plot NMDS in ggplot by Species
plot_sp <- function(scores, centroids){
  # scores - output from get_scores()
  # centroids - output from get_sp_centroids()
  ggplot() + 
    geom_point(data = scores, aes(x = MDS1, y = MDS2, 
                                  color = sp), alpha = .6) + 
    geom_point(data = centroids, aes(x = mean1, y = mean2,
                                     color = sp), size = 3) + 
    geom_errorbar(data = centroids, aes(x = mean1, y = mean2, 
                                        xmin = mean1 - se1, xmax = mean1 + se1,
                                        color = sp), width = 0) + 
    geom_errorbar(data = centroids, aes(x = mean1, y = mean2, 
                                        ymin = mean2 - se2, ymax = mean2 + se2,
                                        color = sp), width = 0) +
    theme_classic(base_size = 14) + 
    scale_color_manual(name = "Species", 
                       values = c("#4381A4", "#48A9A6", "#E4DFDA", "#D4B484", 
                                  "#C1666B", "#896F69","#F09021", "#AA1E57"))
}

# Calculate centroids by Sample Type
get_type_centroids <- function(scores){
  # scores - output from get_scores()
  type.mean1 <- aggregate(scores$MDS1, list(scores$type), mean)
  type.mean2 <- aggregate(scores$MDS2, list(scores$type), mean)
  type.se1 <- aggregate(scores$MDS1, list(scores$type), function(x) sd(x)/sqrt(length(x)))
  type.se2 <- aggregate(scores$MDS2, list(scores$type), function(x) sd(x)/sqrt(length(x)))
  type.dat <- cbind.data.frame(type = type.mean1[1], 
                               mean1 = type.mean1[2], 
                               mean2 = type.mean2[2], 
                               se1 = type.se1[2], 
                               se2 = type.se2[2])
  colnames(type.dat) <- c("type", "mean1", "mean2", "se1", "se2")
  type.dat
}

# Plot NMDS in ggplot by Sample Type
plot_type <- function(scores, centroids){
  # scores - output from get_scores
  # centroids - output from get_type_centroids
  ggplot() + 
    geom_point(data = scores, aes(x = MDS1, y = MDS2, 
                                  color = type), alpha = .6) + 
    geom_point(data = centroids, aes(x = mean1, y = mean2,
                                     color = type), size = 4) + 
    geom_errorbar(data = centroids, aes(x = mean1, y = mean2, 
                                        xmin = mean1 - se1, xmax = mean1 + se1,
                                        color = type), width = 0, size = 1) + 
    geom_errorbar(data = centroids, aes(x = mean1, y = mean2, 
                                        ymin = mean2 - se2, ymax = mean2 + se2,
                                        color = type), width = 0, size = 1) +
    theme_classic(base_size = 14)+ 
    scale_color_manual(name = "Type", 
                       values = c("#4D3371", "#67796F"))
}

##### Estimating Mean Distance -----


# Calculate distance matrix and return long dataframe
get_dist <- function(otu, samps){
  # out - phyloseq@otu_table
  # samps - metadata
  #calculating bray-curtis dissimilarities
  dists <- as.matrix(vegdist(t(otu)))
  dists[lower.tri(dists)] <- NA
  dist.dat <- melt(dists, varnames = c("row", "col"))
  dist.dat <- dist.dat[!is.na(dist.dat$value),]
  names(dist.dat) <- c("samp1", "samp2", "dist")
  # adding sp and type variables
  dist.dat$sp1 <- samps$sp[match(dist.dat$samp1, samps$extract_num)]
  dist.dat$sp2 <- samps$sp[match(dist.dat$samp2, samps$extract_num)]
  dist.dat$type1 <- samps$type[match(dist.dat$samp1, samps$extract_num)]
  dist.dat$type2 <- samps$type[match(dist.dat$samp2, samps$extract_num)]
  dist.dat$same.sp <- ifelse(dist.dat$sp1 == dist.dat$sp2, TRUE, FALSE)
  dist.dat$same.type <- ifelse(dist.dat$type1 == dist.dat$type2, TRUE, FALSE)
  dist.dat
}

# Calculate mean distance per species pair
get_mean_dist <- function(dist.dat){
  #dist.dat - output from get_dist()
  #average distance per species pair
  dist.dat$pair <- interaction(do.call(pmin, dist.dat[4:5]), do.call(pmax, dist.dat[4:5]))
  dist.dat2 <- dist.dat[which(dist.dat$dist != 0), ]
  dist.dat3 <- cbind.data.frame(pair = aggregate(dist.dat2$dist, by = list(dist.dat2$pair), mean)[1], 
                                mean = aggregate(dist.dat2$dist, by = list(dist.dat2$pair), mean)[2], 
                                se = aggregate(dist.dat2$dist, by = list(dist.dat2$pair), function(x) sd(x)/sqrt(length(x)))[2])
  names(dist.dat3) <- c("pair", "mean", "se")
  dist.dat3
}

# Estimate mean distance and return posteriors
get_dist_est <- function(dat){
  # otu - phyloseq object
  stan.fit <- stan(file = "STAN_scripts/meanBC.stan", chains = 4, cores = 4,
                   data = list(D = phyloseq::nsamples(dat),
                               dist = as.matrix(vegdist(t(otu_table(dat)))), 
                               soil_sp = match(sample_data(dat)[['sp']], sp_list$sp),
                               S = length(unique(sample_data(dat)[['sp']]))), 
                   iter = 4000, warmup = 1000)
  as.data.frame(stan.fit)
}

##### Fungi and Competition Models -----

# Make dataframe for long data
uncertain_data <- function(pop_params, dist_params){
  ND_temp <- pivot_longer(data.frame(matrix(pop_params[[3]], nrow= dim(pop_params[[3]])[3], byrow = TRUE)), 
                          names_to = "pair", values_to = "ND", cols = 1:64)
  CR_temp <- pivot_longer(data.frame(matrix(pop_params[[4]], nrow= dim(pop_params[[3]])[3], byrow = TRUE)), 
                          names_to = "pair", values_to = "CR", cols = 1:64)
  DR_temp <- pivot_longer(data.frame(matrix(pop_params[[5]], nrow= dim(pop_params[[3]])[3], byrow = TRUE)), 
                          names_to = "pair", values_to = "DR", cols = 1:64)
  FI_temp <- pivot_longer(data.frame(matrix(pop_params[[6]], nrow= dim(pop_params[[3]])[3], byrow = TRUE)), 
                          names_to = "pair", values_to = "FI", cols = 1:64)
  alpha_temp <- pivot_longer(data.frame(matrix(pop_params[[2]], nrow= dim(pop_params[[3]])[3], byrow = TRUE)), 
                             names_to = "pair", values_to = "alpha", cols = 1:64)
  dist_temp <- pivot_longer(dist_params[,1:64], cols = starts_with("mean"), names_to = "pair", values_to = "dist")
  dat <- cbind.data.frame(ND_temp, CR_temp[,2], DR_temp[,2], FI_temp[,2], alpha_temp[,2], dist_temp[,2], 
                          foc = rep(1:8, 8*nrow(dist_params)), 
                          comp = rep(rep(1:8, each = 8), nrow(dist_params)))
  spp <- cbind.data.frame(num = c(1:8), sp = c("AA", "CA", "CC", "CG", "CP", "NI", "PC", "PN"))
  dat$foc <- spp$sp[match(dat$foc, spp$num)]
  dat$comp <- spp$sp[match(dat$comp, spp$num)]            
  dat$pairs <- interaction(do.call(pmin, dat[8:9]), do.call(pmax, dat[8:9]))
  dat$same <- ifelse(dat$foc == dat$comp, TRUE, FALSE)
  dat$sample <- rep(1:nrow(dist_params), each = 64)
  dat
}

# Fit Asymptotic Models -----

uncertain_curves <- function(params){
  NDintercepts <- NDslopes <- c()
  CRintercepts <- CRslopes <- c()
  DRintercepts <- DRslopes <- lDRintercepts <- lDRslopes <- c()
  FIintercepts <- FIslopes <- lFIintercepts <- lFIslopes <- c()
  for(i in sample(1:nrow(pulls.BH), 5000)){ # indexes for 1000 random samples from posterior
    dat.temp <- params[which(params$sample == i), ]
    dat.temp <- dat.temp[duplicated(dat.temp$pairs), ]
    NDmod.temp <- try(coef(nls(ND ~ a*exp(dist*k), data = dat.temp, 
                               start = list(a = 1, k = 1), algorithm = "port", 
                               lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    lDRmod.temp <- try(coef(nls(log(DR) ~ a*exp(dist*k), data = dat.temp, 
                                start = list(a = 1, k = 1), algorithm = "port", 
                                lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    DRmod.temp <- try(coef(nls(DR - 1 ~ a*exp(dist*k), data = dat.temp, 
                               start = list(a = 1, k = 1), algorithm = "port", 
                               lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    CRmod.temp <- try(coef(nls(CR - 1 ~ a*exp(dist*k), data = dat.temp, 
                               start = list(a = 1, k = 1), algorithm = "port", 
                               lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    lFImod.temp <- try(coef(nls(log(FI) ~ a*exp(dist*k), data = dat.temp, 
                                start = list(a = 1, k = 1), algorithm = "port", 
                                lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    FImod.temp <- try(coef(nls(FI~ a*exp(dist*k), data = dat.temp, 
                               start = list(a = 1, k = 1), algorithm = "port", 
                               lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    NDintercepts <- c(NDintercepts, ifelse(class(NDmod.temp) != "try-error", unname(NDmod.temp[1]), NA))
    NDslopes <- c(NDslopes, ifelse(class(NDmod.temp) != "try-error", unname(NDmod.temp[2]), NA))
    lDRintercepts <- c(lDRintercepts, ifelse(class(lDRmod.temp) != "try-error", unname(lDRmod.temp[1]), NA))
    lDRslopes <- c(lDRslopes, ifelse(class(lDRmod.temp) != "try-error", unname(lDRmod.temp[2]), NA))
    DRintercepts <- c(DRintercepts, ifelse(class(DRmod.temp) != "try-error", unname(DRmod.temp[1]), NA))
    DRslopes <- c(DRslopes, ifelse(class(DRmod.temp) != "try-error", unname(DRmod.temp[2]), NA))
    CRintercepts <- c(CRintercepts, ifelse(class(CRmod.temp) != "try-error", unname(CRmod.temp[1]), NA))
    CRslopes <- c(CRslopes, ifelse(class(CRmod.temp) != "try-error", unname(CRmod.temp[2]), NA))
    lFIintercepts <- c(lFIintercepts, ifelse(class(lFImod.temp) != "try-error", unname(lFImod.temp[1]), NA))
    lFIslopes <- c(lFIslopes, ifelse(class(lFImod.temp) != "try-error", unname(lFImod.temp[2]), NA))
    FIintercepts <- c(FIintercepts, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[1]), NA))
    FIslopes <- c(FIslopes, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[2]), NA))
  }
  cbind.data.frame(NDintercepts, NDslopes, 
                   DRintercepts, DRslopes, 
                   lDRintercepts, lDRslopes,
                   CRintercepts, CRslopes,
                   FIintercepts, FIslopes,
                   lFIintercepts, lFIslopes)
}

# Fit Linear Models -----

uncertain_lines <- function(params){
  ALPHAintercepts <- ALPHAslopes <- INTRAALPHAintercepts <- INTRAALPHAslopes <- INTERALPHAintercepts <- INTERALPHAslopes <- c()
  NDintercepts <- NDslopes <- c()
  CRintercepts <- CRslopes <- c()
  DRintercepts <- DRslopes <- lDRintercepts <- lDRslopes <- c()
  FIintercepts <- FIslopes <- lFIintercepts <- lFIslopes <- c()
  for(i in sample(1:nrow(pulls.BH), 5000)){ # indexes for 1000 random samples from posterior
    dat.temp <- params[which(params$sample == i), ]
    dat.temp2 <- dat.temp[duplicated(dat.temp$pairs), ]
    ALPHAmod.temp <- coef(lm(alpha ~ dist, data = dat.temp))
    INTRAALPHAmod.temp <- coef(lm(alpha ~ dist, data = dat.temp[which(dat.temp$same == TRUE),]))
    INTERALPHAmod.temp <- coef(lm(alpha ~ dist, data = dat.temp[which(dat.temp$same == FALSE),]))
    NDmod.temp <- coef(lm(ND ~ dist, data = dat.temp2))
    lDRmod.temp <- coef(lm(log(DR) ~ dist, data = dat.temp2))
    DRmod.temp <- coef(lm(DR ~ dist, data = dat.temp2))
    CRmod.temp <- coef(lm(CR ~ dist, data = dat.temp2))
    lFImod.temp <- coef(lm(log(FI) ~ dist, data = dat.temp2))
    FImod.temp <- coef(lm(FI ~ dist, data = dat.temp2))
    ALPHAintercepts <- c(ALPHAintercepts, ifelse(class(ALPHAmod.temp) != "try-error", unname(ALPHAmod.temp[1]), NA))
    ALPHAslopes <- c(ALPHAslopes, ifelse(class(ALPHAmod.temp) != "try-error", unname(ALPHAmod.temp[2]), NA))
    INTRAALPHAintercepts <- c(INTRAALPHAintercepts, ifelse(class(INTRAALPHAmod.temp) != "try-error", unname(INTRAALPHAmod.temp[1]), NA))
    INTRAALPHAslopes <- c(INTRAALPHAslopes, ifelse(class(INTRAALPHAmod.temp) != "try-error", unname(INTRAALPHAmod.temp[2]), NA))
    INTERALPHAintercepts <- c(INTERALPHAintercepts, ifelse(class(INTERALPHAmod.temp) != "try-error", unname(INTERALPHAmod.temp[1]), NA))
    INTERALPHAslopes <- c(INTERALPHAslopes, ifelse(class(INTERALPHAmod.temp) != "try-error", unname(INTERALPHAmod.temp[2]), NA))
    NDintercepts <- c(NDintercepts, ifelse(class(NDmod.temp) != "try-error", unname(NDmod.temp[1]), NA))
    NDslopes <- c(NDslopes, ifelse(class(NDmod.temp) != "try-error", unname(NDmod.temp[2]), NA))
    lDRintercepts <- c(lDRintercepts, ifelse(class(lDRmod.temp) != "try-error", unname(lDRmod.temp[1]), NA))
    lDRslopes <- c(lDRslopes, ifelse(class(lDRmod.temp) != "try-error", unname(lDRmod.temp[2]), NA))
    DRintercepts <- c(DRintercepts, ifelse(class(DRmod.temp) != "try-error", unname(DRmod.temp[1]), NA))
    DRslopes <- c(DRslopes, ifelse(class(DRmod.temp) != "try-error", unname(DRmod.temp[2]), NA))
    CRintercepts <- c(CRintercepts, ifelse(class(CRmod.temp) != "try-error", unname(CRmod.temp[1]), NA))
    CRslopes <- c(CRslopes, ifelse(class(CRmod.temp) != "try-error", unname(CRmod.temp[2]), NA))
    lFIintercepts <- c(lFIintercepts, ifelse(class(lFImod.temp) != "try-error", unname(lFImod.temp[1]), NA))
    lFIslopes <- c(lFIslopes, ifelse(class(lFImod.temp) != "try-error", unname(lFImod.temp[2]), NA))
    FIintercepts <- c(FIintercepts, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[1]), NA))
    FIslopes <- c(FIslopes, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[2]), NA))
  }
  cbind.data.frame(ALPHAintercepts, ALPHAslopes,
                   INTERALPHAintercepts, INTERALPHAslopes, 
                   INTRAALPHAintercepts, INTRAALPHAslopes, 
                   NDintercepts, NDslopes, 
                   DRintercepts, DRslopes, 
                   lDRintercepts, lDRslopes,
                   CRintercepts, CRslopes,
                   FIintercepts, FIslopes,
                   lFIintercepts, lFIslopes)
}


# Check histograms of asymptotic model parameter estimates

curvy_histogram <- function(coefs){
  temp_plots <- list()

  temp_plots[["NDslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$NDslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Niche Differences") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$NDslopes)(0), 2), sep = " "))
  
  temp_plots[["NDintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$NDintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Niche Differences") + 
    ylab("Frequency")
  
  temp_plots[["CRslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$CRslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Competitive Ratios") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$CRslopes)(0), 2), sep = " "))
  
  temp_plots[["CRintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$CRintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Competitive Ratios") + 
    ylab("Frequency")
  
  temp_plots[["DRslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$DRslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Demographic Ratios") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$DRslopes)(0), 2), sep = " "))
  
  temp_plots[["DRintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$DRintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Demographic Ratios") + 
    ylab("Frequency")
  
  temp_plots[["lDRslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$lDRslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of log Demographic Ratios") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$lDRslopes)(0), 2), sep = " "))
  
  temp_plots[["lDRintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$lDRintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of log Demographic Ratios") + 
    ylab("Frequency")
  
  temp_plots[["FIslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$FIslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Fitness Inequalities") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$FIslopes)(0), 2), sep = " "))
  
  temp_plots[["FIintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$FIintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Fitness Inequalities") + 
    ylab("Frequency")
  
  temp_plots[["lFIslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$lFIslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of log Fitness Inequalities") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$lFIslopes)(0), 2), sep = " "))
  
  temp_plots[["lFIintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$lFIintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of log Fitness Inequalities") + 
    ylab("Frequency")
  temp_plots
}

# Check histograms of linear model parameter estimates

linear_histogram <- function(coefs){
  temp_plots <- list()
  temp_plots[["ALPHAslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$ALPHAslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Competition") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$ALPHAslopes)(0), 2), sep = " "))
  
  temp_plots[["ALPHAintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$ALPHAintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Competition") + 
    ylab("Frequency")
  
  temp_plots[["INTRAALPHAslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$INTRAALPHAslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Intraspecific Competition") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$INTRAALPHAslopes)(0), 2), sep = " "))
  
  temp_plots[["INTRAALPHAintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$INTRAALPHAintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Intraspecific Competition") + 
    ylab("Frequency")
  
  temp_plots[["INTERALPHAslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$INTERALPHAslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Interspecific Competition") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$INTERALPHAslopes)(0), 2), sep = " "))
  
  temp_plots[["INTERALPHAintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$INTERALPHAintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Interspecific Competition") + 
    ylab("Frequency")
  
  temp_plots[["NDslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$NDslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Niche Differences") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$NDslopes)(0), 2), sep = " "))
  
  temp_plots[["NDintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$NDintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Niche Differences") + 
    ylab("Frequency")
  
  temp_plots[["CRslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$CRslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Competitive Ratios") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$CRslopes)(0), 2), sep = " "))
  
  temp_plots[["CRintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$CRintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Competitive Ratios") + 
    ylab("Frequency")
  
  temp_plots[["DRslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$DRslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Demographic Ratios") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$DRslopes)(0), 2), sep = " "))
  
  temp_plots[["DRintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$DRintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Demographic Ratios") + 
    ylab("Frequency")
  
  temp_plots[["lDRslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$lDRslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of log Demographic Ratios") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$lDRslopes)(0), 2), sep = " "))
  
  temp_plots[["lDRintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$lDRintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of log Demographic Ratios") + 
    ylab("Frequency")
  
  temp_plots[["FIslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$FIslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of Fitness Inequalities") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$FIslopes)(0), 2), sep = " "))
  
  temp_plots[["FIintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$FIintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of Fitness Inequalities") + 
    ylab("Frequency")
  
  temp_plots[["lFIslopes"]] <- ggplot() + 
    geom_histogram(aes(coefs$lFIslopes)) + theme_classic(base_size = 12) + 
    xlab("Growth Rate of log Fitness Inequalities") + 
    ylab("Frequency") + 
    labs(caption = paste("Probability of being < 0:", round(ecdf(coefs$lFIslopes)(0), 2), sep = " "))
  
  temp_plots[["lFIintercepts"]] <- ggplot() + 
    geom_histogram(aes(coefs$lFIintercepts)) + theme_classic(base_size = 12) + 
    xlab("Intercept of log Fitness Inequalities") + 
    ylab("Frequency")
  temp_plots
}

# retreive the medians, 95% CIs, and the % above and below zero for all parameter estimates

just_numbers <- function(coefs){
  temp_list <- list()
  temp_list[["medians"]] <- apply(coefs, 2, median, na.rm = TRUE)
  temp_list[["UCL"]] <- apply(coefs, 2, quantile, 0.975, na.rm = TRUE)
  temp_list[["LCL"]] <- apply(coefs, 2, quantile, 0.025, na.rm = TRUE)
  temp_list[["<0"]] <- apply(coefs, 2, function(x) ecdf(x)(0))
  temp_list[[">0"]] <- apply(coefs, 2, function(x) 1-ecdf(x)(0))
  temp_list[["n"]] <- apply(coefs, 2, function(x) length(x) - sum(is.na(x)))
  temp_list
}

# niche differences plot

curvy_coexistence_points <- function(coefs, points, com){
  ggplot(data = points) + 
    stat_function(size = 1, fun = function(x) median(coefs$NDintercepts, na.rm = TRUE)*exp(x*median(coefs$NDslopes, na.rm = TRUE))) +
    geom_point(aes(x = points[,which(names(points) == com)], y = ND, color = CX), shape = 17, size = 5) +
    scale_color_gradient(low = "#F09021", high = "#48A9A6", name = "Coexistence Probability") + 
    theme_classic(base_size = 12) + #ylim(-1,1) + 
    stat_function(size = .25, fun = function(x) median(coefs$FIintercepts, na.rm = TRUE)*exp(x*median(coefs$FIslopes, na.rm = TRUE))) + 
    # stat_function(size = .25, fun = function(x) (((median(coefs$FIintercepts, na.rm = TRUE)*exp(x*median(coefs$FIslopes, na.rm = TRUE)))) - 1)/
    #                 ((median(coefs$FIintercepts, na.rm = TRUE)*exp(x*median(coefs$FIslopes, na.rm = TRUE))))) + 
    xlab("Fungal Community Dissimilarity") + ylab("Niche Difference")
}

# binomial regression for coexistence probability

CX_model <- function(data){
  CXslopes <- c()
  CXintercepts <- c()
  for(i in sample(1:nrow(pulls.BH), 5000)){ # indexes for 1000 random samples from posterior
    dat.temp <- data[which(data$sample == i), ]
    mod.temp <- glm(CX ~ dist, data = dat.temp, family = "gaussian")
    CXintercepts <- c(CXintercepts, unname(mod.temp[[1]][1]))
    CXslopes <- c(CXslopes, unname(mod.temp[[1]][2]))
  }
  
  cbind.data.frame(CXintercepts, CXslopes, 
                   CXintercepts2 = exp(CXintercepts)/(1+exp(CXintercepts)), 
                   CXslopes2 = exp(CXslopes)/(1+exp(CXslopes)))
}

##### Phylogenetic Analyses -----

uncertain_phylo_data <- function(dist_params, phylo){
  long_dat <- pivot_longer(dist_params[,1:64], cols = starts_with("mean"), names_to = "pair", values_to = "dist")
  long_dat$foc <- substring(long_dat$pair, 11, 11)
  long_dat$comp <- substring(long_dat$pair, 13, 13)
  spp <- cbind.data.frame(num = c(1:8), sp = c("AA", "CA", "CC", "CG", "CP", "NI", "PC", "PN"))
  long_dat$foc2 <- spp$sp[match(long_dat$foc, spp$num)]
  long_dat$comp2 <- spp$sp[match(long_dat$comp, spp$num)]            
  long_dat$pairs <- interaction(do.call(pmin, long_dat[5:6]), do.call(pmax, long_dat[5:6]))
  long_dat$same <- ifelse(long_dat$foc2 == long_dat$comp2, TRUE, FALSE)
  long_dat$p.dist <- phylo$pdist[match(long_dat$pairs, phylo$pair)]
  long_dat$sample <- rep(1:nrow(dist_params), each = 64)
  long_dat
}

uncertain_phylo_lines <- function(params){
  #z-scores
  # params$dist <- scale(params$dist)
  # params$p.dist <- scale(params$p.dist)
  slopes <- c()
  intercepts <- c()
  resids <- list()
  for(i in sample(1:max(params$sample), 5000)){ # indexes for 1000 random samples from posterior
    dat.temp <- params[which(params$sample == i), ]
    mod.temp <- lm(dist ~ p.dist, data = dat.temp)
    intercepts <- c(intercepts, unname(mod.temp[[1]][1]))
    slopes <- c(slopes, unname(mod.temp[[1]][2]))
    resids <- c(resids, resid(mod.temp))
  }
  list(slopes, intercepts, resids)
}

uncertain_phylo_lFI_curves <- function(params, resids){
  FIintercepts <- c()
  FIslopes <- c()
  params <- params[which(params$sample %in% sample(1:max(params$sample), 5000)),]
  for(i in 1:5000){ # indexes for 1000 random samples from posterior
    dat.temp <- params[which(params$sample == i), ]
    res <- as.numeric(resids[i,])
    FImod.temp <- try(coef(nls(log(FI) ~ a*exp(res*k), data = dat.temp, 
                               start = list(a = 1, k = 1), algorithm = "port", 
                               lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    FIintercepts <- c(FIintercepts, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[1]), NA))
    FIslopes <- c(FIslopes, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[2]), NA))
  }
  cbind.data.frame(FIintercepts, FIslopes)
}

uncertain_phylo_FI_curves <- function(params, resids){
  FIintercepts <- c()
  FIslopes <- c()
  params <- params[which(params$sample %in% sample(1:max(params$sample), 5000)),]
  for(i in 1:5000){ # indexes for 1000 random samples from posterior
    dat.temp <- params[which(params$sample == i), ]
    res <- as.numeric(resids[i,])
    FImod.temp <- try(coef(nls(FI - 1 ~ a*exp(res*k), data = dat.temp, 
                               start = list(a = 1, k = 1), algorithm = "port", 
                               lower = c(0, -Inf), upper = c(Inf, Inf))), silent = TRUE)
    FIintercepts <- c(FIintercepts, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[1]), NA))
    FIslopes <- c(FIslopes, ifelse(class(FImod.temp) != "try-error", unname(FImod.temp[2]), NA))
  }
  cbind.data.frame(FIintercepts, FIslopes)
}


uncertain_phylo_ND_curves <- function(params, resids){
  NDintercepts <- c()
  NDslopes <- c()
  params <- params[which(params$sample %in% sample(1:max(params$sample), 5000)),]
  for(i in 1:5000){ # indexes for 5000 random samples from posterior
    dat.temp <- params[which(params$sample == i), ]
    res <- as.numeric(resids[i,])
    NDmod.temp <- try(coef(nls(ND ~ 1 + a*exp(res*k), data = dat.temp, 
                               start = list(a = -1, k = -1), algorithm = "port", 
                               lower = c(-Inf, -Inf), upper = c(0, Inf))), silent = TRUE)
    NDintercepts <- c(NDintercepts, ifelse(class(NDmod.temp) != "try-error", unname(NDmod.temp[1]), NA))
    NDslopes <- c(NDslopes, ifelse(class(NDmod.temp) != "try-error", unname(NDmod.temp[2]), NA))
  }
  cbind.data.frame(NDintercepts, NDslopes)
}


