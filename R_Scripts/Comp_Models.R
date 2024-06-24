### Fitting Competition Models

set.seed(6)

library(rstan)
library(abind)
library(dplyr)
library(loo)
library(MCMCvis)
library(tidybayes)
library(tidyverse)

dat2 <- rbind.data.frame(dat[which(dat$comp != "A"), ], 
                         data.frame(id = rep(dat$id[which(dat$comp == "A")], 8), 
                                    foc = rep(dat$foc[which(dat$comp == "A")], 8), 
                                    comp = c(rep("AA", length(dat$id[which(dat$comp == "A")])), 
                                             rep("CA", length(dat$id[which(dat$comp == "A")])), 
                                             rep("CC", length(dat$id[which(dat$comp == "A")])), 
                                             rep("CG", length(dat$id[which(dat$comp == "A")])),
                                             rep("CP", length(dat$id[which(dat$comp == "A")])), 
                                             rep("NI", length(dat$id[which(dat$comp == "A")])), 
                                             rep("PC", length(dat$id[which(dat$comp == "A")])), 
                                             rep("PN", length(dat$id[which(dat$comp == "A")]))), 
                                    fec = rep(dat$fec[which(dat$comp == "A")], 8), 
                                    count = rep(dat$count[which(dat$comp == "A")], 8)))


# a sp list to convert sp names to numbers
sp_list <- cbind.data.frame(num = 1:length(unique(dat2$foc)), sp = sort(unique(dat2$foc)))

# saveRDS(sp_list, "RDS_files/sp_list.RDS")

foc_sp <- match(dat2$foc, sp_list$sp)
comp_sp <- match(dat2$comp, sp_list$sp)
n_plants <- length(foc_sp)
n_species <- length(unique(dat2$foc))
comp_count <- dat2$count
fec <- dat2$fec

stan_dat <- list(sp = foc_sp,
                 comp = comp_sp, 
                 count = comp_count, 
                 N = n_plants, 
                 S = n_species, 
                 fec = fec)

# fitting four population models:
# null: F = λ
# Linear Lotka-Volterra: F = λ-αij*Nj
# Beverton-Holt: F = λ/(1 + αij*Nj)
# Ricker: F = λe^(-αij*Nj)

# un-comment these lines if you do not have the RDS files
# 
# stan.fitnull <- stan(file = "STAN_scripts/null_stan_model.stan", data = stan_dat, chains = 4, cores = 4,
#                      iter = 4000, warmup = 1000, pars = c("lambda", "log_lik"))
# stan.fitLV <- stan(file = "STAN_scripts/LV_stan_model.stan", data = stan_dat, chains = 4, cores = 4,
#                    iter = 4000, warmup = 1000, pars = c("lambda", "alpha", "ND", "CR", "DR", "FI", "log_lik")) # some convergence issues here
stan.fitBH <- stan(file = "STAN_scripts/BH_stan_model.stan", data = stan_dat, chains = 4, cores = 4,
                 iter = 4000, warmup = 1000, pars = c("lambda", "alpha", "sigma", "ND", "CR", "DR", "FI", "log_lik"))
# stan.fitRick <- stan(file = "STAN_scripts/Ricker_stan_model.stan", data = stan_dat, chains = 4, cores = 4,
#                      iter = 4000, warmup = 1000, pars = c("lambda", "alpha", "ND", "CR", "DR", "FI", "log_lik"))


# saveRDS(stan.fitnull, file = "RDS_files/NULLfit.RDS")
# saveRDS(stan.fitLV, file = "RDS_files/LVfit.RDS")
saveRDS(stan.fitBH, file = "RDS_files/BHfit.RDS")
# saveRDS(stan.fitRick, file = "RDS_files/Rickfit.RDS")

# reading in the RDS files that contain the rstan objects generated + saved in previous lines

stan.fitnull <- readRDS("RDS_files/NULLfit.RDS")
stan.fitLV <- readRDS("RDS_files/LVfit.RDS")
stan.fitRick <- readRDS("RDS_files/Rickfit.RDS")
stan.fitBH <- readRDS("RDS_files/BHfit.RDS")

# saving some pdfs of the trace plots
# MCMCtrace(stan.fitnull, params = "lambda", filename = "null_trace.pdf")
# MCMCtrace(stan.fitLV, params = c("lambda", "alpha"), filename = "LV_trace.pdf")
# MCMCtrace(stan.fitBH, params = c("lambda", "alpha"), filename = "BH_trace.pdf")
# MCMCtrace(stan.fitRick, params = c("lambda", "alpha"), filename = "Rick_trace.pdf")

# Model Comparison via LOO
comps <- loo_compare(list(loo(extract_log_lik(stan.fitnull)), 
                          loo(extract_log_lik(stan.fitLV)), 
                          loo(extract_log_lik(stan.fitBH)), 
                          loo(extract_log_lik(stan.fitRick))))

comps # BH fits best, so I'll proceed with those estimates
