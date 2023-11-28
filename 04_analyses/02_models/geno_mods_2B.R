### This script will explore plant performance relationships from the GCR
### task 2A experiment to rank the genotypes as more or less "interesting" for
### future phases of the experiment. 

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

biomass <- read.csv("data/biomassdata.csv")
census <- read.csv("data/GCR_2A_CensusData.csv")
licor<- read.csv("data/LiCorData/licordataclean.csv")

### Clean the data -------------------------------------------------------------
licor2 <- licor %>% group_by(pot) %>%
  summarise(meanPSN = mean(Photo),
            meanTPN = mean(Trmmol),
            leafarea = mean(leafarea)) %>%
  mutate(PSN = meanPSN/leafarea,
         TPN = meanTPN/leafarea)

dat <- left_join(census, biomass, by = "pot")
dat0<- left_join(dat, licor2, by = "pot")
#View(dat)

dat1<- dat0 %>% filter(CensusNo == 1)  %>% #list of OG plants from exp
  mutate(stat = if_else(massA > 0, 1, 0),
         stat = coalesce(stat,0),
         Soil = str_extract(pot, "\\w"),
         Geno = str_extract(pot, "\\w{2,}"),
         Level = substr(pot, 3,3),
         massT = massA + massB) %>%
  #remove pots that don't have plants in them
  filter(Soil != "C" &
           Soil != "P" &
           Soil != "N" &
           Soil != "S") #ok I should have just done Soil == "L"...

outs <- c("515", "639", "A20", "D16", "R24", "SWE", "212", "105")
#genotypes that have not survived well

dat2 <- dat1 %>% filter(!Geno %in% outs)

### First question- which genos are most likely to survive? --------------------

dat2$Geno <- as.factor(dat2$Geno)
dat2$Level <- as.factor(dat2$Level)
contrasts(dat2$Geno) <- contr.sum(9)
contrasts(dat2$Level) <- contr.sdif(3)

mod1 <- glmer(stat ~ Geno + Level + log(HtPrevCensus) + (1|table),
               family = binomial,
               data = dat2)
summary(mod1)  

saveRDS(mod1, file = "04_analyses/02_models/Output/SurvGeno1.RDS")
#Levels: CRI G15 G74 INA K19 MOC TAS VIR YON
#MOC and TAS are slightly more likely to survive overall
#INA has slightly lower probability of survival
# surivaly is sequentially more likely with higher soil moisture

mod2 <- glmer( stat ~ Geno*Level + (1|table),
               family = binomial,
               data = dat2)
summary(mod2)
#well that's trash, but not surprising. 

#based on mod1 and my notes the following are candidate genos based on survival
#(filter #1); TAS, MOC, CRI, YON, G15, VIR, K19

### Second question- are there any interesting biomass patterns among the 
### remaining genos?

outs2<- c("G74", "INA")
dat3 <- dat2 %>% filter(!Geno %in% outs2)

## total biomass:
dat3$Geno <- as.factor(dat3$Geno)
contrasts(dat3$Geno) <- contr.sum(7)
mod3<- lmer(log(massT) ~ Geno + Level + (1|table),
            data = dat3)
summary(mod3)

mod4<- lmer(log(massA/massB) ~ Geno + Level + (1|table),
            data = dat3)
summary(mod4)

mod5<- lmer(log(massB) ~ Geno + Level + (1|table),
            data = dat3)
summary(mod5)

#TAS, MOC, and CRI are definite top 3 genos. MOC and TAS are good at surviving, 
# but TAS has smaller than average roots. CRI has larger roots than average. 

### third question- interesting physiological patterns?

mod6<- lmer(meanPSN ~ Geno + Level + (1|table),
            data = dat3)
summary(mod6) #nothing

mod7<- lmer(meanTPN ~ Geno + Level + (1|table),
            data = dat3)
summary(mod7) #nothing
