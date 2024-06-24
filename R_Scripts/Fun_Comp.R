### Fungi and Competition Analysis

set.seed(6)

library(tidyverse)
library(gridExtra)

##### Bringing data in -----

# population parameter estimates

pulls.BH <- as.data.frame(readRDS("RDS_files/BHfit.RDS"))

# mean fungal distance estimates
pulls.tot <- readRDS("RDS_files/estDist_tot.RDS")
pulls.patho <- readRDS("RDS_files/estDist_patho.RDS")
pulls.sapro <- readRDS("RDS_files/estDist_sapro.RDS")
pulls.symbio <- readRDS("RDS_files/estDist_symbio.RDS")

# extracting pop params and calculating coexistence metrics

lambdas <- pulls.BH[,1:8]
alphas <- pulls.BH[,9:72]

alphas2 <- array(data = NA, dim = c(8, 8, nrow(alphas)))
for(i in 1:nrow(alphas)){
  alphas2[,,i] <- as.numeric((alphas[i,]))
}

ND <- CR1 <- CR2 <- CR <- DR1 <- DR2 <- DR <- index <- FI <- array(data = NA, dim = c(8, 8, nrow(alphas)))
for(f in 1:8){
  for(c in 1:8){
    ND[f,c,] <- 1/sqrt((alphas2[f,c,]*alphas2[c,f,])/(alphas2[f,f,]*alphas2[c,c,]))
    CR1[f,c,] <- sqrt((alphas2[c,c,]*alphas2[c,f,])/(alphas2[f,f,]*alphas2[f,c,]))
    CR2[f,c,] <- sqrt((alphas2[f,f,]*alphas2[f,c,])/(alphas2[c,c,]*alphas2[c,f,]))
    DR1[f,c,] <- lambdas[,f]/lambdas[,c] #might need to add - 1's
    DR2[f,c,] <- lambdas[,c]/lambdas[,f]
    index <- ifelse(CR1*DR1 > CR2*DR2, 1, 2)
    FI <- ifelse(CR1*DR1 > CR2*DR2, CR1*DR1, CR2*DR2)
    CR <- ifelse(index == 1, CR1, CR2)
    DR <- ifelse(index == 1, DR1, DR2)
  }
}

params <- list(lambdas, alphas2, ND, CR, DR, FI)

# saveRDS(params, file = "RDS_files/pop_params")

med.alpha <- apply(alphas, 2, median)
med.alpha2 <- apply(alphas2, c(1,2), median)
med.ND <- apply(ND, c(1,2), median)
med.FI <- apply(FI, c(1,2), median)
med.CR <- apply(CR, c(1,2), median)
med.DR <- apply(DR, c(1,2), median)

CX_prob <- apply(ifelse(FI < ND, 1, 0), c(1,2), function(x) sum(x)/length(x))
EX_prob <- apply(ifelse(FI > ND & FI > 1/(1-ND), 1, 0), c(1,2), function(x) sum(x)/length(x))
PE_prob <- apply(ifelse(FI < 1 - ND & FI > 1/(1-ND), 1, 0), c(1,2), function(x) sum(x)/length(x))

# extracting median fungal dissimilarities

dists <- cbind.data.frame(rbind(pulls.tot, pulls.patho, pulls.sapro, pulls.symbio), 
                          com = rep(c("tot", "patho", "sapro", "symbio"), each = nrow(pulls.tot)))

dist_points <- list()
for(i in c("tot", "patho", "sapro", "symbio")){
  mat.med <- matrix(apply(dists[which(dists$com == i), grepl("mean", names(dists))], 2, median), nrow = 8, byrow = FALSE)
  mat.med2 <- matrix(0, ncol = ncol(mat.med), nrow = nrow(mat.med))
  mat.med2[lower.tri(mat.med2, diag = TRUE)] <- mat.med[lower.tri(mat.med2, diag = TRUE)] 
  dist_points[[i]][["median"]] <- c(mat.med2 + t(mat.med2) - diag(diag(mat.med2))) # ensuring that the matrix is symmetric
  
  mat.low <- matrix(apply(dists[which(dists$com == i), grepl("mean", names(dists))], 2, quantile, 0.05), nrow = 8, byrow = FALSE)
  mat.low2 <- matrix(0, ncol = ncol(mat.low), nrow = nrow(mat.low))
  mat.low2[lower.tri(mat.low2, diag = TRUE)] <- mat.low[lower.tri(mat.low2, diag = TRUE)]
  dist_points[[i]][["lower"]] <- matrix(c(mat.low2 + t(mat.low2) - diag(diag(mat.low2))), nrow = 8, ncol = 8, )
  
  mat.up <- matrix(apply(dists[which(dists$com == i), grepl("mean", names(dists))], 2, quantile, 0.95), nrow = 8, byrow = FALSE)
  mat.up2 <- matrix(0, ncol = ncol(mat.up), nrow = nrow(mat.up))
  mat.up2[lower.tri(mat.up2, diag = TRUE)] <- mat.up[lower.tri(mat.up2, diag = TRUE)]
  dist_points[[i]][["upper"]] <- c(mat.up2 + t(mat.up2) - diag(diag(mat.up2)))
}

# compiling sp pair dataset

sp_list <- read_rds("RDS_files/sp_list.RDS")
sp_pair_dat <- cbind.data.frame(foc = which(upper.tri(med.ND), arr.ind = TRUE)[,1], 
                                comp = which(upper.tri(med.ND), arr.ind = TRUE)[,2], 
                                ND = med.ND[upper.tri(med.ND)], FI = med.FI[upper.tri(med.FI)], 
                                CR = med.CR[upper.tri(med.CR)], DR = med.DR[upper.tri(med.DR)], 
                                tot = dist_points$tot$median[which(upper.tri(med.ND))], 
                                patho = dist_points$patho$median[which(upper.tri(med.ND))], 
                                sapro = dist_points$sapro$median[which(upper.tri(med.ND))], 
                                symbio = dist_points$symbio$median[which(upper.tri(med.ND))], 
                                CX = CX_prob[which(upper.tri(CX_prob))], 
                                EX = EX_prob[which(upper.tri(EX_prob))], 
                                PE = PE_prob[which(upper.tri(PE_prob))])

sp_pair_dat$foc <- sp_list$sp[match(sp_pair_dat$foc, sp_list$num)]
sp_pair_dat$comp <- sp_list$sp[match(sp_pair_dat$comp, sp_list$num)]
sp_pair_dat$pair <- interaction(do.call(pmin, sp_pair_dat[1:2]), do.call(pmax, sp_pair_dat[1:2]))

# saveRDS(sp_pair_dat, file = "RDS_files/sp_pair_dat.RDS")

# compiling competition dataset

comp_dat <- cbind.data.frame(foc = substr(names(med.alpha), 7, 7), comp = substr(names(med.alpha), 9, 9), 
                             alpha = as.numeric(med.alpha), tot = dist_points$tot$median, 
                             patho = dist_points$patho$median, sapro = dist_points$sapro$median, 
                             symbio = dist_points$symbio$median)

comp_dat$foc <- sp_list$sp[match(comp_dat$foc, sp_list$num)]
comp_dat$comp <- sp_list$sp[match(comp_dat$comp, sp_list$num)]
comp_dat$type <- ifelse(comp_dat$foc == comp_dat$comp, "intra", "inter")

### Competition ~ Fungi -----

# compiling dataset containing all pulls from posterior distributions
comp.tot <- uncertain_data(params, pulls.tot)
comp.patho <- uncertain_data(params, pulls.patho)
comp.sapro <- uncertain_data(params, pulls.sapro)
comp.symbio <- uncertain_data(params, pulls.symbio)

# asymptotic curves for ND, CR, DR, and FI

# this takes awhile... only run if you don't have the RDSs
curves.tot <- uncertain_curves(comp.tot)
curves.patho <- uncertain_curves(comp.patho)
curves.sapro <- uncertain_curves(comp.sapro)
curves.symbio <- uncertain_curves(comp.symbio)
write_rds(curves.tot, "RDS_files/curves.tot.RDS")
write_rds(curves.patho, "RDS_files/curves.patho.RDS")
write_rds(curves.sapro, "RDS_files/curves.sapro.RDS")
write_rds(curves.symbio, "RDS_files/curves.symbio.RDS")
curves.tot <-readRDS("RDS_files/curves.tot.RDS")
curves.patho <-readRDS("RDS_files/curves.patho.RDS")
curves.sapro <-readRDS("RDS_files/curves.sapro.RDS")
curves.symbio <-readRDS("RDS_files/curves.tot.RDS")

# linear models for intra-alpha, inter-alpha, ND, CR, DR, and FI

lines.tot <- uncertain_lines(comp.tot)
lines.patho <- uncertain_lines(comp.patho)
lines.sapro <- uncertain_lines(comp.sapro)
lines.symbio <- uncertain_lines(comp.symbio)

# for looking at the distribution of parameter estimates
histogram_list <- list()
j <- 0
for(i in list(curves.tot, curves.patho, curves.sapro, curves.symbio)){
  j <- j + 1
  histogram_list[["curves"]][[j]] <- curvy_histogram(i)
}
j <- 0
for(i in list(lines.tot, lines.patho, lines.sapro, lines.symbio)){
  j <- j + 1
  histogram_list[["lines"]][[j]] <- linear_histogram(i)
}

# for checking the descriptive stats for the parameter estimates
numbers_list <- list()
j <- 0
for(i in list(curves.tot, curves.patho, curves.sapro, curves.symbio)){
  j <- j + 1
  numbers_list[["curves"]][[j]] <- just_numbers(i)
}
j <- 0
for(i in list(lines.tot, lines.patho, lines.sapro, lines.symbio)){
  j <- j + 1
  numbers_list[["lines"]][[j]] <- just_numbers(i)
}

# plotting

curvy_coexistence_points(curves.tot, sp_pair_dat, "tot") + 
  theme(legend.position = "none") # Figure 3
curvy_coexistence_points(curves.patho, sp_pair_dat, "patho") + 
  theme(legend.position = "none")
curvy_coexistence_points(curves.sapro, sp_pair_dat, "sapro") + 
  theme(legend.position = "none")
curvy_coexistence_points(curves.symbio, sp_pair_dat, "symbio") + 
  theme(legend.position = "none")

### Coexistence Outcome ~ Fungi -----

names(sp_pair_dat)[length(names(sp_pair_dat))] <- "pairs"
cx_tot <- merge(comp.tot, sp_pair_dat[,c("CX", "pairs")], by = "pairs")
cx_patho <- merge(comp.patho, sp_pair_dat[,c("CX", "pairs")], by = "pairs")
cx_sapro <- merge(comp.sapro, sp_pair_dat[,c("CX", "pairs")], by = "pairs")
cx_symbio <- merge(comp.symbio, sp_pair_dat[,c("CX", "pairs")], by = "pairs")

CXmod.tot <- CX_model(cx_tot)
CXmod.patho <- CX_model(cx_patho)
CXmod.sapro <- CX_model(cx_sapro)
CXmod.symbio <- CX_model(cx_symbio)

coefs <- cbind.data.frame(CXmod.tot$CXslopes, CXmod.tot$CXintercepts)
lines_dat<- apply(coefs, 1, function(z) {
  stat_function(fun = function(x)exp(z[1]*x + z[2])/(1 + exp(z[1]*x + z[2])), colour="darkgrey")})

# Figure 4
ggplot(data = sp_pair_dat, aes(x = tot, y = CX)) +
  lines_dat[sample(1:length(lines_dat), 50)] + 
  geom_point(size = 2) + 
  theme_classic(base_size = 15) + 
  stat_function(linewidth = 1.25, 
                fun = function(x) exp(median(CXmod.tot$CXslopes)*x + median(CXmod.tot$CXintercepts))/(1 + exp(median(CXmod.tot$CXslopes)*x + median(CXmod.tot$CXintercepts)))) + 
  xlab("Fungal Community Dissimilarity") + ylab("Probability of Coexistence")



