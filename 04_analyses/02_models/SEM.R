rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(MASS) #for contr.sdif
library(broom.mixed) #model diagnostics
library(performance) #model diagnostics
library(ggeffects) # for direct easy model prediction plotting
library(piecewiseSEM) # for SEM

dat <- read.csv("data/ModData/FullCensusTimeSeries.csv")
RGR <- readRDS("04_analyses/02_models/Output/RGR.RDS")
Biomass <- readRDS("04_analyses/02_models/Output/Biomass.RDS")

#RGR, Biomass, and Survival all use different datasets so you can't do an SEM
#with them without some pretty coding. 

#Instead, what are the geochemical relationships of interest here?
# I'm starting to think that this might not be super necessary... brain is tired though...