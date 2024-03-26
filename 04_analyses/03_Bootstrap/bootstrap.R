### Bootstrap the model output

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(ggplot2) # plot results
library(patchwork) #joins plots nicely
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

### Read in Models -------------------------------------------------------------
#plant performance
TotBM <- readRDS("04_analyses/02_models/Output/Biomass.RDS")
RMF <- readRDS("04_analyses/02_models/Output/RMF.RDS")
Mort <- readRDS("04_analyses/02_models/Output/Mortality.RDS")
RGR <- readRDS("04_analyses/02_models/Output/RGR.RDS")

### Bootstrap the Model --------------------------------------------------------
##Performance -------------------------------------------

#[1]
TotBM.boot <- bootMer(TotBM, nsim = 1000, FUN = function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         Level =c("1", "2", "3"))
  predict(., newdata = preddat, re.form = ~0)
},
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)

#[2]
RMF.boot <- bootMer(RMF, nsim = 1000, FUN = function(.) {
  preddat <- expand.grid(totBM = 0.089,
                         Soil = c("L", "S"),
                         Level =c("1", "2", "3"))
  predict(., newdata = preddat, re.form = ~0)
},
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)

#[3]
RGR.boot <- bootMer(RGR, nsim = 1000, FUN = function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         Level =c("1", "2", "3"),
                         HtPrevCensus = 280)
  predict(., newdata = preddat, re.form = ~0)
},
                    parallel="snow", ncpus = detectCores(), 
                    cl=cl)
#[4]
#this isn't working for some reason. Why? 
#Mort.boot <- bootMer(Mort, nsim = 1000, FUN = function(.) {
 # preddat <- expand.grid(Soil = c("L", "S"),
  #                       Level =c("1", "2", "3"),
  #                       HtPrevCensus = 280,
  #                       censusint = 14)
#  predict(., newdata = preddat, re.form = ~0)
#},
 #                   parallel="snow", ncpus = detectCores(), 
  #                  cl=cl)
#stop the clusters once bootstrapping is finished-------
stopCluster(cl = cl) 

### Save the Output ------------------------------------------------------------
##Performance ---------------------------------------
save(TotBM.boot, file = "04_analyses/03_Bootstrap/output/TotBM.boot")
save(RMF.boot, file = "04_analyses/03_Bootstrap/output/RMF.boot")
save(RGR.boot, file = "04_analyses/03_Bootstrap/output/RGR.boot")
#save(Mort.boot, file = "04_analyses/03_Bootstrap/output/Mort.boot")
