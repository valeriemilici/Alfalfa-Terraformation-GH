#bartlett test, levene's test

### Initialize workspace -------------------------------------------------------
rm(list =ls())
library(lme4)
library(tidyverse) #clean data manipulation
library(car)
library(parallel)
library(stats)

pH <- readRDS("04_analyses/02_models/Output/pHSoil.RDS") 
ECmod <- readRDS("04_analyses/02_models/Output/ECSoil.RDS")

# Generate Predicted Data (remove effect of Level) 
pred.fun <- function(pH) {
  model.matrix(~Soil,
               data = expand.grid(Soil =c("L","C", "N", "S"))) %*%
    fixef(pH)[c(1:4)]
}

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))
##EC plot and Bartlett ---------------------------------------------------------
ECBoot <- bootMer(ECmod, FUN = pred.fun,
                nsim = 1000, parallel = "snow", ncpus = detectCores(),cl = cl)
#Prepare Data 
ECdraws <- data.frame(ECBoot$t)

Soil <- c("L", "C", "N", "S")

names(ECdraws) <- Soil

#for plotting
test <- pivot_longer(ECdraws, everything())

#Bartlett Test 
bartlett.test(list(exp(ECdraws$S), exp(ECdraws$L)))
#Weird. Without the effect of soil moisture everything is different from each
#other, except L and S are similar to each other... gotta make a figure and see 
#what's up

#plot of model predictions
ggplot(test, aes(x = name, y = exp(value))) + geom_boxplot() + theme_bw()


##pH plot and Bartlett ---------------------------------------------------------
pHBoot <- bootMer(pH, FUN = pred.fun,
                nsim = 1000, parallel = "snow", ncpus = detectCores(),cl = cl)
#Prepare Data 
pHdraws <- data.frame(pHBoot$t)

Soil <- c("L", "C", "N", "S")

names(pHdraws) <- Soil

#for plotting
pHdat <- pivot_longer(pHdraws, everything())

#Bartlett Test 
bartlett.test(list((pHdraws$C), (pHdraws$N)))
#Everything has a different distribution from everything else. L and N look
#pretty similar to me though

leveneTest(value ~ name, pHdat)

#plot of model predictions
ggplot(pHdat, aes(x = name, y = (value))) + geom_boxplot() + theme_bw()
