
##### MODEL FITTING CHECK #####

library(bayesplot)
library(rstanarm)
library(scales)

source("R_scripts/Data_Compilation.R")

mod <- readRDS("RDS_files/BHfit.RDS")
mod.df <- as.data.frame(mod)

dat$foc <- match(dat$foc, unique(dat$foc))
dat$comp <- match(dat$comp, unique(dat$comp)[-1])
predictions <- matrix(data = NA, nrow = nrow(mod.df), ncol = nrow(dat))
for(i in 1:nrow(mod.df)){
  lambdas <- mod.df[i,1:8]
  alphas <- matrix(data = mod.df[i,9:(8+64)], 
                   ncol = 8, nrow = 8, byrow = FALSE)
  sigma <- mod.df[i,(8+64+1):(8+64+1+8)]
  
  for(j in 1:nrow(dat)){
    mu <- lambdas[dat$foc[j]]/(1 + 
                                 ifelse(is.na(dat$comp[j]), 0,
                                        dat$count[j]*alphas[[dat$foc[j],dat$comp[j]]]))
    pred <- rnorm(1, mean = as.numeric(mu), sd = as.numeric(sigma[dat$foc[j]]))
    predictions[i,j] <- pred
  }
}

saveRDS(predictions, file = "RDS_files/pp_check")

pp_check_df <- cbind.data.frame(med_pred = apply(predictions, 2, median), 
                                low_pred = apply(predictions, 2, quantile, 0.025), 
                                up_pred = apply(predictions, 2, quantile, 0.975), 
                                real = dat$fec, 
                                sp = dat$foc,
                                par = 1:ncol(predictions))

sum(pp_check_df$real > pp_check_df$low_pred & pp_check_df$real < pp_check_df$up_pred)/
  nrow(pp_check_df)

ggplot(data = pp_check_df[pp_check_df$sp == 1,], aes(x = par, y = med_pred, 
                               ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  labs(title = "Acmispon americanus") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 2,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(title = "Clarkia amoena") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 3,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  labs(title = "Calindrinia ciliata") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 4,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  labs(title = "Collomia grandiflora") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 5,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  labs(title = "Clarkia purpurea") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 6,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  labs(title = "Navarretia intertexta") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 7,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(title = "Plectritis congesta") +
  theme(plot.title = element_text(face = "italic"))

ggplot(data = pp_check_df[pp_check_df$sp == 8,], aes(x = par, y = med_pred, 
                                                     ymin = low_pred, ymax = up_pred)) + 
  geom_point() + geom_errorbar(width = 0) + 
  geom_point(aes(x = par, y = real), color = "red") + 
  theme_classic(base_size = 15) + ylab("Predictions") + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(title = "Plagiobothrys nothofulvus") +
  theme(plot.title = element_text(face = "italic"))

# mean difference between predictions and real values?

meds <- lows <- ups <- c()
for(i in 1:ncol(predictions)){
  x <- predictions[,i] - dat$fec[i]
  meds <- c(meds, median(x))
  lows <- c(lows, quantile(x, 0.025))
  ups <- c(ups, quantile(x, 0.975))
}
mean 
ggplot(data = cbind.data.frame(par = 1:ncol(predictions), 
                               med_dif = meds, 
                               low_dif = lows, 
                               up_dif = ups), 
       aes(x = par, y = med_dif, ymin = low_dif, ymax = up_dif)) + 
  geom_point() + geom_errorbar()


