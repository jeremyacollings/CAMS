### Compilation of Competition Data

set.seed(6)

library(tidyverse)
library(readxl)


### Aggregating all of the fecundity data ----

AAseeds <- read_excel("SeedsData.xlsx", sheet = "ACMAME")
AAstructures <- read_excel("StructuresData.xlsx", sheet = "ACMAME")

MeanSeeds <- mean(AAseeds$N_Seeds)

ACMAME_dat <- data.frame(id = AAstructures$ID, 
                         foc = AAstructures$Foc, 
                         comp = AAstructures$Comp, 
                         dens = AAstructures$Dens, 
                         fec = AAstructures$N_fruits*MeanSeeds)
ACMAME_dat[ACMAME_dat == "-"] <- "A"

ACMAME_dat$dens <- factor(ACMAME_dat$dens, levels = c("A", "L", "M", "H"))

ACMAME_dat$dens <- ifelse(ACMAME_dat$dens == "A", 0, 
                          ifelse(ACMAME_dat$dens == "L", 1, 
                                 ifelse(ACMAME_dat$dens == "M", 2, 3)))


CCseeds <- read_excel("SeedsData.xlsx", sheet = "CALCIL")
CCflowers <- read_excel("FlowerData.xlsx", sheet = "CALCIL")
CCstructures <- read_excel("StructuresData.xlsx", sheet = "CALCIL")

CCseeds[CCseeds == "M2"] <- "M"

mean(CCseeds$N_Seeds)
median(CCseeds$N_Seeds)

MeanSeeds <- mean(CCseeds$N_Seeds)

mean(CCflowers$N_flow)
median(CCflowers$N_flow)

MeanFlowers <- mean(CCflowers$N_flow)

CALCIL_dat <- data.frame(id = CCstructures$ID, 
                         foc = CCstructures$Foc, 
                         comp = CCstructures$Comp, 
                         dens = CCstructures$Dens, 
                         fec = CCstructures$N_Stems*MeanSeeds*MeanFlowers)
CALCIL_dat[CALCIL_dat == "-"] <- "A"

CALCIL_dat$dens <- factor(CALCIL_dat$dens, levels = c("A", "L", "M", "H"))


CALCIL_dat$dens <- ifelse(CALCIL_dat$dens == "A", 0, 
                          ifelse(CALCIL_dat$dens == "L", 1, 
                                 ifelse(CALCIL_dat$dens == "M", 2, 3)))


CAseeds <- read_excel("SeedsData.xlsx", sheet = "CLAAMO")
CAstructures <- read_excel("StructuresData.xlsx", sheet = "CLAAMO")


mean(CAseeds$N_Seeds)
median(CAseeds$N_Seeds)
MeanSeeds <- mean(CAseeds$N_Seeds)

CLAAMO_dat <- data.frame(id = CAstructures$ID, 
                         foc = CAstructures$Foc, 
                         comp = CAstructures$Comp, 
                         dens = CAstructures$Dens, 
                         fec = CAstructures$N_Fruit*MeanSeeds)
CLAAMO_dat[CLAAMO_dat == "-"] <- "A"


CLAAMO_dat$dens <- ifelse(CLAAMO_dat$dens == "A", 0, 
                          ifelse(CLAAMO_dat$dens == "L", 1, 
                                 ifelse(CLAAMO_dat$dens == "M", 2, 3)))

CPseeds <- read_excel("SeedsData.xlsx", sheet = "CLAPUR")
CPstructures <- read_excel("StructuresData.xlsx", sheet = "CLAPUR")

mean(CPseeds$N_Seeds)
median(CPseeds$N_Seeds)

MeanSeeds <- mean(CPseeds$N_Seeds)

CLAPUR_dat <- data.frame(id = CPstructures$ID, 
                         foc = CPstructures$Foc, 
                         comp = CPstructures$Comp, 
                         dens = CPstructures$Dens, 
                         fec = CPstructures$N_Fruit*MeanSeeds)
CLAPUR_dat[CLAPUR_dat == "-"] <- "A"


CLAPUR_dat$dens <- ifelse(CLAPUR_dat$dens == "A", 0, 
                          ifelse(CLAPUR_dat$dens == "L", 1, 
                                 ifelse(CLAPUR_dat$dens == "M", 2, 3)))

CGseeds <- read_excel("SeedsData.xlsx", sheet = "COLLOM")
CGflowers <- read_excel("FlowerData.xlsx", sheet = "COLLOM")
CGstructures <- read_excel("StructuresData.xlsx", sheet = "COLLOM")

CGseeds$N_Seeds[CGseeds$N_Seeds == "f"] <- 4

MeanSeeds <- median(as.numeric(CGseeds$N_Seeds))


MeanFlowers <- median(CGflowers$N_flow)

COLLOM_dat <- data.frame(id = CGstructures$ID, 
                         foc = CGstructures$Foc, 
                         comp = CGstructures$Comp, 
                         dens = CGstructures$Dens, 
                         fec = CGstructures$N_Flow_Head*MeanSeeds*MeanFlowers)
COLLOM_dat[COLLOM_dat == "-"] <- "A"

COLLOM_dat$dens <- factor(COLLOM_dat$dens, levels = c("A", "L", "M", "H"))

NIseeds <- read_excel("SeedsData.xlsx", sheet = "NAVINT")
NIflowers <- read_excel("FlowerData.xlsx", sheet = "NAVINT")
NIstructures <- read_excel("StructuresData.xlsx", sheet = "NAVINT")

MeanSeeds <- mean(NIseeds$N_Seeds)


mean(NIflowers$N_flow)
median(NIflowers$N_flow)

MeanFlowers <- median(NIflowers$N_flow)

NAVINT_dat <- data.frame(id = NIstructures$ID, 
                         foc = NIstructures$Foc, 
                         comp = NIstructures$Comp, 
                         dens = NIstructures$Dens, 
                         fec = NIstructures$N_Heads*MeanSeeds*MeanFlowers)
NAVINT_dat[NAVINT_dat == "-"] <- "A"

NAVINT_dat$dens <- factor(NAVINT_dat$dens, levels = c("A", "L", "M", "H"))

PNseeds <- read_excel("SeedsData.xlsx", sheet = "PLANOT")
PNflowers <- read_excel("FlowerData.xlsx", sheet = "PLANOT")
PNstructures <- read_excel("StructuresData.xlsx", sheet = "PLANOT")

MeanSeeds <- mean(PNseeds$N_Seeds)


PNflowers$Size[PNflowers$Size == "-"] <- "LG"

MeanFlowersSmall <- mean(PNflowers$N_flow[PNflowers$Size == "SM"])
MeanFlowersLarge <- mean(PNflowers$N_flow[PNflowers$Size == "LG"])

PLANOT_dat <- data.frame(id = PNstructures$ID, 
                         foc = PNstructures$Foc, 
                         comp = PNstructures$Comp, 
                         dens = PNstructures$Dens, 
                         fec = (PNstructures$N_LG_Stems*MeanSeeds*MeanFlowersLarge +
                                  PNstructures$N_SM_Stems*MeanSeeds*MeanFlowersSmall))

PLANOT_dat$dens <- factor(PLANOT_dat$dens, levels = c("A", "L", "M", "H"))

PCflowers <- read_excel("PC_FlowerData.xlsx")
PCstructures <- read_excel("PC_StructuresData.xlsx")


mean(PCflowers$count)
MedianFlowers <- median(PCflowers$count)


PLECON_dat <- data.frame(id = PCstructures$ID, 
                         foc = PCstructures$Foc, 
                         comp = PCstructures$Comp, 
                         dens = PCstructures$Dens, 
                         fec = PCstructures$Count*MedianFlowers)

fecund_dat <- rbind.data.frame(ACMAME_dat, CALCIL_dat, CLAAMO_dat, 
                               CLAPUR_dat, COLLOM_dat, NAVINT_dat, 
                               PLANOT_dat, PLECON_dat)

fecund_dat[nrow(fecund_dat) + 1,] <- c("AA_CA_H1", "AA", 
                                       "CA", 3, 
                                       sum(fecund_dat[12:13, 5]))
fecund_dat <- fecund_dat[-c(12, 13), ]
fecund_dat[which(fecund_dat$comp == "CL"), ] <- c("CP_CA_L1", "CP", 
                                                  "CA", 1, 
                                                  fecund_dat[150, 5])

fecund_dat[which(fecund_dat$foc == "Ca"), ] <- c("CA_A1_2", "CA", 
                                                 "A", 0, 
                                                 fecund_dat[116, 5])
fecund_dat[which(fecund_dat$foc == "CN"), ] <- c("CP_PN_M1", "CP", 
                                                 "PN", 2, 
                                                 fecund_dat[158, 5])

fecund_dat$comp[which(fecund_dat$comp == "-")] <- "A"

### Competitor data ----

backgrounds <- read_excel("backgrounds.xlsx")

backgrounds2 <- backgrounds[backgrounds$id %in% fecund_dat$id, ]

`%!in%` = Negate(`%in%`)
fecund_dat$id[which(fecund_dat$id %!in% backgrounds2$id)]
backgrounds$id[which(backgrounds$id %!in% fecund_dat$id)]

dat <- merge(fecund_dat[, -4], backgrounds2, by = "id")
dat <- dat[, c(1:4, 8)]
dat$counts[which(dat$comp == "A" & 
                                   is.na(dat$counts))] <- 0
names(dat) <- c("id", "foc", "comp", "fec", "count")
dat[which(dat$comp == "CF"), ] <- c("AA_CG_L1", "AA", "CG", 
                                    dat[24, 4], dat[24, 5])
dat$fec <- as.numeric(dat$fec)
dat$count <- as.numeric(dat$count)

View(dat)
