### This is an analysis on the 2A alfalfa data that asks...
### is there a difference in  survival between seedlings growing in
### live vs. sterile & by moisture/genos?

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(broom.mixed) #model diagnostics
library(performance) #model diagnostics
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

dat <- read.csv("data/performance.csv") 

### clean the data--------------------------------------------------------------
outs <- c("515", "639", "A20", "D16", "R24", "SWE", "212", "105", "G74")
#genotypes that have not survived well

dat1 <- dat %>% filter(!Geno %in% outs)
dat1$HtPrevCensus[dat1$CensusNo == 1] <- 0 #set Census 1 prev ht to 0 not NA
dat1 <- dat1 %>% 
  group_by(pot) %>% 
  arrange(CensusNo) %>%
  filter(row_number() == n())

dat2 <- dat %>% filter(!Geno %in% outs,
                       CensusNo >= 2 & HtPrevCensus > 0)
dat3 <- dat %>%  filter(!Geno %in% outs) %>% group_by(pot) %>% 
  arrange(CensusNo) %>%
  filter(row_number() == n())
#thinned data set, row is removed after a plant is marked as dead the first time

### set the contrasts ----------------------------------------------------------

dat3$Geno <- as.factor(dat3$Geno)
dat3$Geno <- factor(dat3$Geno,
                        levels = c("CRI", "INA", "K19", "MOC",
                                   "TAS", "VIR", "YON", "G15"))
dat2$Level <- as.factor(dat2$Level)
contrasts(dat3$Geno) <- contr.sum(8) #intercept represents ave across genos
dat3$Level <- as.factor(dat3$Level)
contrasts(dat3$Level) <- contr.sdif(3)
contrasts(dat2$Level) <- contr.helmert(3) #intercept is ave across soil water

### model 1 survival by soil and moisture --------------------------------------
dat2$Soil <- factor(dat2$Soil, levels = c("S", "L"))
Modsurv <- glmer(status ~ Soil + Level  + scale(HtPrevCensus) +
                  (1|pot)+ (1|CensusNo) + (1|Geno) + (1|table),
                 data = dat2,
                 family = binomial(link = "cloglog"),
                 glmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 10000)),
                 offset = scale(log(censusint)))

check_model(Modsurv)

summary(Modsurv)

# Soil and level are very correlated and can't interact in the model. 
# when the model is additive we see that plants are more likely to survive in
# sterile soil rather than live soil (p = 0.003) and that plants are more likely
# to survive as soil moisture increases (p = 0.006 & p = 0.013). Obviously
# plants that are taller are more likely to survive. And as we saw in RGR mods
# plants in sterile soil and at higher moisture grew faster == were taller ==
# are more likely to survive to the next census. 

### model 2 soil and water survival begin and end only -------------------------
dat1$Soil <- factor(dat1$Soil, levels = c("L", "S"))

Modsurv2 <- glmer(status ~ Soil + Level + (1|Geno) + (1|table),
                 data = dat3,
                 family = binomial(),
                 glmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 10000)))
#all pots need a census 6 reading of 0 or 1, check.
check_model(Modsurv)

summary(Modsurv2)

### Genotype Specific ----------------------------------------------------------
dat1$Soil <- factor(dat1$Soil, levels = c("L", "S"))

survGeno <- glmer(status ~ Level + Geno + (1|table),
                  data = filter(dat3, Soil == "L"),
                  family = binomial(),
                  glmerControl(optimizer = "bobyqa",
                               optCtrl = list(maxfun = 10000)))
summary(survGeno)

