### This is an analysis on the 2A alfalfa data that asks...
### is there a difference in final performance metrics between seedlings growing
###in live vs. sterile & by moisture/genos?

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(MASS)
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
PSNdat <- dat %>% filter(PSN >= 0) %>%
  group_by(pot) %>%
  arrange(CensusNo) %>%
  filter(row_number() == n())

PSNdat1 <- dat1 %>% filter(PSN >= 0) %>%
  group_by(pot) %>%
  arrange(CensusNo) %>%
  filter(row_number() == n())

#set contrasts------------------------------------------------------------------

PSNdat$Level <- as.factor(PSNdat$Level)
contrasts(PSNdat$Level) <- contr.sdif(3)

PSNdat1$Level <- as.factor(PSNdat1$Level)
contrasts(PSNdat1$Level) <- contr.sdif(3)

PSNdat1$Geno <- factor(PSNdat1$Geno,
                      levels = c("CRI", "INA", "K19", "MOC",
                                 "TAS", "VIR", "YON", "G15"))
contrasts(PSNdat1$Geno) <- contr.sum(8)

# [1] Photosynthesis Soil * Level ----------------------------------------------
PSNSoil <- lmer(PSN ~ Soil*Level + (1|Geno) + (1|table),
                data = filter(PSNdat, PSN < 0.07))

check_model(PSNSoil) #decent

plot(PSNSoil)
summary(PSNSoil) #same story as always. No interactions.

# [2] Photosynthesis Genotype --------------------------------------------------

PSNGeno <- lmer(PSN ~ Soil + Level*Geno + (1|table),
                data = filter(PSNdat1, PSN < 0.07),
                REML = F)

check_model(PSNGeno) #decent

plot(PSNGeno)
summary(PSNGeno) #VIR experienced a reduction in photosynthesis when soil 
#moisture increased. 
saveRDS(PSNGeno, file = "04_analyses/02_models/Output/PSNGeno.RDS")
# [3] Transpiration Soil * Level -----------------------------------------------
TPNSoil <- lmer(TPN ~ Soil*Level + (1|Geno) + (1|table),
                data = filter(PSNdat, TPN < 0.16))

check_model(TPNSoil)
plot(TPNSoil)
summary(TPNSoil) #only thing to report is that sterile plants reduced
#transpiration at the highest soil moisture levels. So the least transpiration
#when soil moisture was greatest and there was plenty of water. Curious. 

# [4] Transpiration Genos ------------------------------------------------------

TPNGeno <- lmer(TPN ~ Soil + Level*Geno + (1|table),
                data = filter(PSNdat1, TPN < 0.11),
                REML = F)

check_model(TPNGeno)
plot(TPNGeno)
summary(TPNGeno) #Yon has above average transpiration

saveRDS(TPNGeno, file = "04_analyses/02_models/Output/TPNGeno.RDS")
