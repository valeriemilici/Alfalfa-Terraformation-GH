### Linear Models with PCA output

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

#PCA output data
perfdat <- read.csv("data/perfcoords1.csv")
chemdat <- read.csv("data/discoords.csv")
#Treatment/location data
meta <- read.csv("data/metadata.csv")

### clean the data--------------------------------------------------------------
outs <- c("515", "639", "A20", "D16", "R24", "SWE", "212", "105", "G74")
#genotypes that have not survived well

#note- I am leaving perfdat and chemdat as separate dataframes instead of 
#joining them, because perfdat is longer, and I don't want to lose data.
#I could join them and have many NAs, but I don't want to do that.

#add pot and treatment information to the data
perfdat1 <- left_join(perfdat, meta, by = "pot") 
#for chemdat you first need to rename the column
chemdat1<- rename(chemdat, pot = "dat5.pot")
chemdat1 <- left_join(chemdat1, meta, by = "pot")
#remove the random x.x and x.y columns
perfdat1 <- perfdat1 %>% select(!c(X.x, X.y))
chemdat1<- chemdat1 %>% select(!c(X.x, X.y))

#remove the failed genotypes
perfdat2 <- perfdat1 %>% filter(!Geno %in% outs)
chemdat2 <- chemdat1 %>% filter(!Geno %in% outs)

#fix some NAs
chemdat2$Geno <- chemdat2$Geno %>% replace_na("CON")
chemdat2$Soil <- chemdat2$Soil %>% replace_na("C")
chemdat2$Level <- chemdat2$Level %>% replace_na(2)
### The models -----------------------------------------------------------------

## Set the Contrasts ------------------------
perfdat1$Level <- as.factor(perfdat1$Level)
perfdat2$Level <- as.factor(perfdat2$Level)
chemdat1$Level <- as.factor(chemdat1$Level)
chemdat2$Level <- as.factor(chemdat2$Level)
contrasts(perfdat1$Level) <- contr.sdif(3)
contrasts(perfdat2$Level) <- contr.sdif(3)
contrasts(chemdat1$Level) <- contr.sdif(3)
contrasts(chemdat2$Level) <- contr.sdif(3)

perfdat2$Geno <- factor(perfdat2$Geno,
                      levels = c("CRI", "INA", "K19", "MOC",
                                 "TAS", "VIR", "YON", "G15"))

chemdat2$Geno <- factor(chemdat2$Geno,
                        levels = c("CRI", "INA", "K19", "MOC",
                                   "TAS", "VIR", "YON", "CON", "G15"))
chemdat2$Soil <- factor(chemdat2$Soil,
                        levels = c("L", "S", "N", "C"))
contrasts(perfdat2$Geno) <- contr.sum(8)
contrasts(chemdat2$Geno) <- contr.sum(9)

## Performance PCA models -------------------

#Effect of genotype, etc. on PC1 -------
mod1 <- lmer(Dim.1 ~ Soil + Geno*Level + (1|Location),
             data = perfdat2)

check_model(mod1)
#no collinearity problems, overall model looks decent.
#maybe under predicts the high values

summary(mod1)
#Sterile soil is associated with higher values of PC1, which makes sense bc
#PC1 is all about plant size. No genotype deviates within the "size cluster"
#of PC1, except for K19, which is generally "small". There is no interaction
#between genotypes and soil moisture. In general plants are larger between 
#level 1 and level 2, but there is not difference between level 2 and 3. 
saveRDS(mod1, file = "04_analyses/02_models/Output/perfPCAD1.RDS")

#effect of genotype on PC2 ------
mod2 <- lmer(Dim.2 ~ Soil + Geno*Level + (1|Location),
             data = perfdat2)

check_model(mod2)
#there is one really high value that makes things weird. Otherwise looks pretty
#good.

summary(mod2)
#no difference in gas exchange between sterile and live soil plants. Interesting.
#VIR has less gas exchange than the average geno. We see that at level 2, CRI
#increases its gas exchange. INA seems to have above average gas exchange at level 3
# and TAS, interestingly, has worse gas-exchange at level 3. 

saveRDS(mod2, file = "04_analyses/02_models/Output/perfPCAD2.RDS")
## Geochem PCA models -------------------------

#Effect of genotype on PC1 ----
mod3 <- lmer(Dim.1 ~ Soil + Geno + Level + (1|Location),
             data = chemdat2)

check_model(mod3)
#this model looks considerably better without the interaction term

summary(mod3)
#Plants in Live soil have higher PC1 values (Chlorine, Sulfate, pH, EC) than 
#plants in sterile soil. N/C soils don't differ from each other, so it's a 
#plants + microbes story. They must be together. MOC is associated with higher
#PC1 values. Also, PC1 is lower at the highest soil moisture level. Of this, 
# the model is quite certain. 

saveRDS(mod3, file = "04_analyses/02_models/Output/chemPCAD1.RDS")
#Effect of genotype on PC2 ------
mod4 <- lmer(Dim.2 ~ Soil + Geno + Level + (1|Location),
             data = chemdat2)

check_model(mod4)
#this model looks mediocre. The distribution of the data are slightly non-normal.
#I should look into narrower distributions to try to fix this. But there's 
#prob nothing in this model anyway so let's check it out!

summary(mod4)
#PC2 is NO3 and F. Only MOC shows any pattern, and that's lower values of PC2. 
#there are no other patterns here. Perhaps a better model fit would help. 
saveRDS(mod4, file = "04_analyses/02_models/Output/chemPCAD2.RDS")
