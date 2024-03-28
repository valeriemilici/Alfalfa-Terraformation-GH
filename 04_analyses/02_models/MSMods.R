## The models specifically used in the manuscript. If modified from the main
## Analysis report/markdown document

rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(performance) #model diagnostics
library(stats) #bootstrap
library(parallel) #bootstrap in parallel
library(ggeffects) #for quick beautiful plotting

dat <- read.csv("data/ModData/AllPerformanceGeochem.csv")

dat <- filter(dat, Soil != "P")

dat$Geno <- dat$Geno %>% replace_na("C")

dat$Level <- as.factor(dat$Level)

### Trends in Porewater Geochemistry

dat1 <- dat %>% mutate(soil_CNP = case_when( Soil == "L" ~ "P",
                                             Soil == "S" ~ "P",
                                             .default = Soil))
#Sodium - model to create figure 6A
mod6A <- lmer(log(mu_Na) ~ soil_CNP * Level + scale(Dmass_tot) +
                 (1|Geno) + (1|dis_table) + (1|table) + (1|batch_cat),
               data = dat1)

saveRDS(mod6A, "04_analyses/02_models/Output/mod6A.RDS")

#Lithium- model to create figure 6B

mod6B <- lmer(mu_Li ~ soil_CNP * Level + scale(Dmass_tot) +
                 (1|Geno) + (1|dis_table) + (1|table) + (1|batch_cat),
               data = filter(dat1, mu_Li < 1.5))
saveRDS(mod6B, "04_analyses/02_models/Output/mod6B.RDS")
#Magnesium - model to create figure 6C

mod6C <- lmer(log(mu_Mg) ~ Soil * Level + scale(Dmass_tot) +
                 (1|Geno) + (1|dis_table) + (1|table) + (1|batch_cat),
               data = dat)

saveRDS(mod6C, "04_analyses/02_models/Output/mod6C.RDS")


