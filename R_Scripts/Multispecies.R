
##### MULTISPECIES ANALYSES FOR CAMS #####

struc_dat <- rbind.data.frame(readRDS("CAMS_struc_out1.RDS"), 
                                  readRDS("CAMS_struc_out2.RDS"), 
                                  readRDS("CAMS_struc_out3.RDS"), 
                                  readRDS("CAMS_struc_out4.RDS"))

struc_dat <- struc_dat[,-ncol(struc_dat)]
struc_dat$rich <- nchar(gsub("_", "", struc_dat$combos))
struc_dat$iter <- rep(1:4000, each = 247)

percents <- richness <- c()
for(i in 1:4000){
  x <- struc_dat[struc_dat$iter == i,]
  percents <- c(percents, 
                unname(tapply(x$feasabilities, x$rich, function(x) sum(x)/length(x))))
  richness <- c(richness, 
                as.numeric(names(tapply(x$feasabilities, x$rich, function(x) sum(x)/length(x)))))
}

tapply(percents, richness, median)
tapply(percents, richness, quantile, 0.025)
tapply(percents, richness, quantile, 0.975)

ggplot(data = cbind.data.frame(rich = names(tapply(percents, richness, median)), 
                               med = tapply(percents, richness, median), 
                               low = tapply(percents, richness, quantile, 0.025), 
                               up = tapply(percents, richness, quantile, 0.975)), 
       aes(x = rich, y = med, ymin = low, ymax = up)) + 
  geom_point(size = 2) + geom_errorbar(width = 0) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  theme_classic(base_size = 15) + xlab("Species Richness") + 
  ylab("Percent Combinations Feasible")

CPD_df <- cbind.data.frame(med = tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                  struc_dat$combos[which(struc_dat$rich != 2)], median), 
                           low = tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                        struc_dat$combos[which(struc_dat$rich != 2)], quantile, 0.025), 
                           up = tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                       struc_dat$combos[which(struc_dat$rich != 2)], quantile, 0.975), 
                           combo = names(tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                                struc_dat$combos[which(struc_dat$rich != 2)], quantile, 0.025)))
CPD_df$rich <- nchar(gsub("_", "", CPD_df$combo))

sum(CPD_df$low < 0 & CPD_df$up > 0)/nrow(CPD_df)

ggplot(data = CPD_df, aes(x = fct_reorder(combo, rich), y = med, 
                          ymin = low, ymax = up, color = as.factor(rich))) + 
  geom_point() + geom_errorbar(width = 0, alpha = .5) + 
  theme_classic(base_size = 15) + ylab("Community Pair Differential") + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  xlab("Species Combination") + 
  scale_y_continuous(breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)) + 
  scale_color_manual(name = "Species Richness", 
                     values = c("#4381A4", "#F09021", "#E4DFDA", 
                                "#D4B484", "#48A9A6", "#AA1E57"))

ggplot(data = cbind.data.frame(rich = names(tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                                   struc_dat$rich[which(struc_dat$rich != 2)], median)), 
                               med = tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                            struc_dat$rich[which(struc_dat$rich != 2)], median), 
                               low = tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                            struc_dat$rich[which(struc_dat$rich != 2)], quantile, 0.025), 
                               up = tapply(struc_dat$differentials[which(struc_dat$rich != 2)], 
                                           struc_dat$rich[which(struc_dat$rich != 2)], quantile, 0.975)), 
       aes(x = rich, y = med, ymin = low, ymax = up)) + 
  geom_point(size = 2) + geom_errorbar(width = 0) + 
  theme_classic(base_size = 15) + xlab("Species Richness") + 
  ylab("Community-Pair Differential") + 
  geom_hline(yintercept = 0, linetype = "dashed")

ggplot(data = struc_dat, 
       aes(x = as.factor(rich), y = coex_rates)) + 
  geom_violin()


ggplot(data = struc_dat, 
       aes(x = combos, y = coex_rates)) + 
  geom_violin()

ggplot(data = struc_dat, 
       aes(x = combos, y = omegas)) + 
  geom_violin() + facet_wrap(~ rich)

