### This is an analysis on the 2A alfalfa data that asks...
### is there a difference in final performance metrics between seedlings growing
###in live vs. sterile & by moisture/genos?

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(MASS) #for contr.sdif
library(broom.mixed) #model diagnostics
library(performance) #model diagnostics
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

dat <- read.csv("data/performance.csv") 

### clean the data--------------------------------------------------------------
outs <- c("515", "639", "A20", "D16", "R24", "SWE", "212", "105", "G74")
#genotypes that have not survived well

dat1 <- dat %>% filter(!Geno %in% outs) 

#Remove NAs and simplify data to one observation per plant
BMdat <- dat %>% filter(massA > 0) %>%
  group_by(pot) %>%
  arrange(CensusNo) %>%
  filter(row_number() == n()) %>%
  mutate(totBM = massA + massB,
         BMrat = massB/massA)

BMdat1 <- dat1 %>% filter(massA > 0) %>%
  group_by(pot) %>%
  arrange(CensusNo) %>%
  filter(row_number() == n()) %>%
  mutate(totBM = massA + massB)

#set contrasts------------------------------------------------------------------

BMdat$Level <- as.factor(BMdat$Level)
contrasts(BMdat$Level) <- contr.sdif(3)

BMdat1$Level <- as.factor(BMdat1$Level)
contrasts(BMdat1$Level) <- contr.sdif(3)

BMdat1$Geno <- factor(BMdat1$Geno,
                      levels = c("CRI", "INA", "K19", "MOC",
                                 "TAS", "VIR", "YON", "G15"))
contrasts(BMdat1$Geno) <- contr.sum(8)

#[1] Aboveground + failed genos ------------------------------------------------
BMSoilA <- lmer(log(massA) ~ Soil * Level + (1|Geno) + (1|table),
                data = BMdat)

check_model(BMSoilA)
#indicates some problem with collinearity when soil and level interact.
#indicates row 162 and 64 may be overly influential

ggplot(data.frame(lev=hatvalues(BMSoilA),pearson=residuals(BMSoilA,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw() 

levId <- which(hatvalues(BMSoilA) >= .12)

BMSoilA1 <- lmer(log(massA) ~ Soil * Level + (1|Geno) + (1|table),
                data = BMdat[-c(as.vector(levId)), ])

check_model(BMSoilA1)

summary(BMSoilA1) #model A1 looks pretty good, removing the points does not 
#affect the inference. It is better to include the failed genotypes for this
#model

#[2] Belowground + failed genos------------------------------------------------

#When including roots, remove obs that had no root sample taken. Will increase
#their root mass relative to others and skew data.

BMSoilB<- lmer(log(massB) ~ log(totBM) + Soil*Level + (1|Geno) + (1|table),
               data = filter(BMdat, rootsample == "Y"))

check_model(BMSoilB) #the model looks pretty good. Consider a t-distribution
#there appears to be some fat tails. 

ggplot(data.frame(lev=hatvalues(BMSoilB),pearson=residuals(BMSoilB,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw() 
#looks pretty good and clustered

summary(BMSoilB) #plants have larger roots in sterile soil & overall roots get 
#larger as it get wetter. 

#[3] Biomass Ratio ------------------------------------------------------------

BMSoilR<- lmer(log(massA/massB) ~ log(totBM) + Soil*Level + (1|Geno) + (1|table),
                   data = filter(BMdat, rootsample == "Y"))
check_model(BMSoilR) #again this model looks ok. There may be some influential
#points (look into interpreting Cook's distance) but otherwise this is a solid 
#model.

summary(BMSoilR) 
# Overall plants invest more in belowground biomass than aboveground biomass and
# this is more pronounced in plants in sterile soil (p < 0.001). However, as
# soil moisture increases, there is a little jump for sterile plants, where they
# invest slightly more in aboveground allocation (but plants are still mainly 
# belowground invested) (p = 0.04). Could this be evidence that something was
# attacking/consuming the roots of the plants in live soil?
saveRDS(BMSoilR, file = "04_analyses/02_models/Output/BMSoilR.RDS")

#[4] Total biomass ------------------------------------------------------------
BMdat$Soil <- factor(BMdat$Soil, levels = c("L", "S"))

#BMdat$totBM.c <- BMdat$totBM - mean(BMdat$totBM)

BMSoilT<- lmer(log(totBM) ~ Soil*Level + (1|Geno) + (1|table),
                   data = filter(BMdat, rootsample == "Y"),
               REML = F)

check_model(BMSoilT) #Once more a decent fit. Same problems as above. Biomass
# allocation model has the best fit. 

summary(BMSoilT) #it's the tale we know and love. Sterile plants are bigger. all
#plants like it wetter.
saveRDS(BMSoilT, file = "04_analyses/02_models/Output/BMSoilT.RDS")

BMdat2 <- BMdat %>% group_by(Soil, Level) %>%
  summarise(aveBM = mean(totBM))
#[5] genotype-specific ratio --------------------------------------------------
BMGenoR<- lmer(log(massA/massB) ~ log(totBM) + Soil + Level*Geno + (1|table),
               data = filter(BMdat1, rootsample == "Y"))

check_model(BMGenoR) 
#this looks good! Again Cook's distance could be fixed maybe

summary(BMGenoR)
# MOC and TAS invest more in aboveground biomass allocation compared to the
# average genotype. YON invests more in belowground allocation compared to the
# average genotype. Soil moisture does not affect these patterns. 
saveRDS(BMGenoR, file = "04_analyses/02_models/Output/BMGenoR.RDS")
#[6] genotype-specific total biomass ------------------------------------------

BMGenoT<- lmer(log(totBM) ~ Soil + Level*Geno + (1|table),
               data = filter(BMdat1, rootsample == "Y"))

 check_model(BMGenoT) #not to shabby!

summary(BMGenoT) #K19 is slightly smaller than average (p = 0.11). MOC gets
#larger as soil moisture increases (0.10). VIR gets smaller as soil moisture
# increases (p = 0.048). 
saveRDS(BMGenoT, file = "04_analyses/02_models/Output/BMGenoT.RDS")


#[7] genotype-specific root biomass ------------------------------------------

BMGenoR<- lmer(log(massB) ~ log(totBM) +Soil + Level*Geno + (1|table),
               data = filter(BMdat1, rootsample == "Y"))

check_model(BMGenoT) #not to shabby!

summary(BMGenoR) #MOC has slightly smaller belowground biomass (p < 0.001).
#But it increases as soil moisture increases (p = 0.04) and overcomes that previous 
#decrease.
#TAS has also slightly smaller belowground biomass allocation (p = 0.004). 

#[8] genotype-specific root biomass ------------------------------------------

BMGenoS<- lmer(log(massA) ~ log(totBM) +Soil + Level*Geno + (1|table),
               data = filter(BMdat1, rootsample == "Y"))

check_model(BMGenoT) #not to shabby!

summary(BMGenoS) #MOC and TAS have bigger shoots than the average geno (p <0.001
# & p = 0.033), and YON has smaller shoots than the average geno (p = 0.01). 
#moisture does not affect this relationship. 



