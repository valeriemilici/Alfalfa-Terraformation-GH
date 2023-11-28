### Discharge Models -----
### Water Geochemistry (anions, pH, and ec)

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(MASS) #cotr.sdif
library(broom.mixed) #model diagnostics
library(performance) #model diagnostics
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

dat <- read.csv("data/2A_porewater_geochemistry.csv") 
perf <- read.csv("data/performance.csv")
pcadat<- read.csv("data/discoords.csv")
dis <- read.csv("data/discharge.csv")
location <- read.csv("data/GCR_2A_CensusData.csv")
ic <- read.csv("data/NPOC.csv")

### Data manipulation-----------------------------------------------------------

##For PCAs-----------------------
#rename column so it matches other data
colnames(pcadat)[4] <- "pot"
#merge geochemistry pca output with performance data
perf1 <- left_join(perf, pcadat, by = "pot") %>%
  filter(!is.na(Dim.1), !is.na(Dim.2), CensusNo == 6) %>%
  mutate(biorat = massA/massB)
##For all models -------------------
#extract location data for random effect
location2 <- location %>% filter(CensusNo == 6) %>%
  mutate(block = paste(table,rack, sep = "."))
location2 <- location2[,c(3,7,12)] #Why isn't select working?!

dat1 <- left_join(dat, location2, by = "pot")
dat1 <- left_join(dat1, ic, by = "pot")

#Set Contrasts -----
# Geochem Data ------------------
dat1$Level <- as.factor(dat1$Level)
contrasts(dat1$Level) <- contr.sdif(3)

dat1$Geno <- dat1$Geno %>% replace_na("control")
dat1$Geno <- factor(dat1$Geno)

dat1$soil <- factor(dat1$soil,
                   levels = c("L", "S", "B"))

contrasts(dat1$Geno) <- contr.sum(12)

dat1$TC_batch <- as.factor(dat1$TC_batch)
contrasts(dat1$TC_batch) <- contr.sum(3)
# Performance Data ---------------------
perf1$Level <- as.factor(perf1$Level)
contrasts(perf1$Level) <- contr.sdif(3)

perf1$Geno <- as.factor(perf1$Geno)
contrasts(perf1$Geno) <- contr.sum(11)

## Genotype-Specific Patterns --------
outs <- c("515", "639", "A20", "D16", "R24", "SWE", "212", "105", "G74")
#genotypes that have not survived well

dat2 <- dat1 %>% filter(!Geno %in% outs) 

perf2 <- perf1 %>% filter(!Geno %in% outs)

dat2$Geno <- factor(dat2$Geno,
                    levels = c("CRI", "K19", "MOC",
                               "TAS", "VIR", "G15", "YON"))
contrasts(dat2$Geno) <- contr.sum(7)

perf2$Geno <- factor(perf2$Geno,
                     levels = c("CRI", "K19", "MOC",
                                "TAS", "VIR", "G15", "YON"))
contrasts(perf2$Geno) <- contr.sum(7)

## Soil levels ---------------------------
dat1$Soil <- relevel(factor(dat1$Soil), 
                     ref = "C")
dat1$Geno <- relevel(factor(dat1$Geno),
                     ref = "MOC")
#pH ------------------------------------------------------------------------
pHSoil <- lmer(pH ~ Soil * Level + scale(Dmass_tot) + (1|Geno) + (1|block),
                data = dat1)


check_model(pHSoil)
#indicates row 162 and 64 may be overly influential

summary(pHSoil)
#Soil treatment has basically no effect on pH.
#Intermediate soil moisture is associated with an increase in pH
# (2-1 = 0.585, 0.143, p < 0.001), while high soil moisture is associated with
# a decrease in pH (3-2 = -0.462, 0.109, p < 0.001. As you might expect, 
# as discharge volume increase, pH decreases (-0.291, 0.046, p < 0.001)
# One difference is that at soil moisture level #2, plants in sterile soil have
# lower pore-water pH than plants in live soil (-0.411, 0.169, 0.012). 

saveRDS(pHSoil, "04_analyses/02_models/Output/pHSoil.RDS")

#pHGeno <- lm(pH ~ Soil + Geno * Level,
#               data = dat2)
#check_model(pHGeno)
#summary(pHGeno)

#this is pretty messy. move on. (try to scale the data...)

#EC ------------------------------------------------------------------------
ECSoil <- lmer(log(EC) ~ Soil * Level + scale(Dmass_tot) + (1|Geno) + (1|block),
               data = dat1)
#logging looks a little better for EC.
check_model(ECSoil)

summary(ECSoil)
#nothing affects EC other than discharge volume

#save model output
saveRDS(ECSoil, file = "04_analyses/02_models/Output/ECSoil.RDS")

#EC x biomass 
dat4<- left_join(dat1, perf, by = "pot")
dat4 <- dat4 %>% group_by(pot) %>% arrange(CensusNo) %>%
  filter(row_number() == n()) %>% 
  mutate(totBM = massA + massB,
         biorat = massA/massB)

ECbiomass <- lmer(log(EC) ~ log(totBM)+ Level.x + (1|Geno.x) + (1|table),
                  data = dat4)
summary(ECbiomass) #biomass has absolutely no effect on EC (total, above, and below)

ECGeno <- lmer(log(EC) ~ Soil + Geno * Level + scale(Dmass_tot) + (1|block),
             data = dat2)
check_model(ECGeno)
summary(ECGeno)
#There is no variation in EC
saveRDS(ECGeno, file = "04_analyses/02_models/Output/ECGeno.RDS")
#Fluorine? -----------------------------------------------------------------

FSoil <- lmer(log(Fluorine) ~ Soil * Level + (1|Geno) + (1|block),
               data = dat1)

check_model(FSoil)

summary(FSoil)
#this looks terrible. Also the data have two peaks so I'm not sure that this
#data meets the assumptions I am trying to make

##If necessary, try analyzing fluorine with a GAM

# Nitrate ----------------------------------------------------------------------

NO3Soil <- lmer(Nitrate ~ Soil*Level + scale(Dmass_tot) + (1|Geno) + (1|block),
                data = dat1)

check_model(NO3Soil)

summary(NO3Soil) 
#There may be more nitrate in the abiotic (C) compared to the full biotic (L) 
# treatments (18.43, 9.15, 0.098). N pots also have lower nitrate than C pots
# (-17.687, 8.01, 0.031). Thus, microbes, not plants lead to a decrease in
# pore-water nitrate. For live plants, nitrate also increases between the dry and
# intermediate treatments, but not between the intermediate and high 
# (2-1 = 19.02, 9.14, 0.041). At higher soil moisture levels, sterile plants (S) have
# lower nitrate than live plants (-23.107, 10.28, 0.027). 

saveRDS(NO3Soil, file = "04_analyses/02_models/Output/NO3Soil.RDS")

#genotype specific
NO3Geno <- lmer(log(Nitrate) ~ Soil + Geno*Level + scale(Dmass_tot) + (1|block),
              data = dat2)
check_model(NO3Geno)
summary(NO3Geno)
#TAS has below-average nitrate (-0.909, 0.436, 0.044), it gets even lower when
# soil moisture increases (2-1 = -2.13, 0.77, 0.009).
# Meanwhile, MOC has above-average nitrate when soil moisture increases
# marginal (1.226, 0.707, 0.092)
saveRDS(NO3Geno, file = "04_analyses/02_models/Output/NO3Geno.RDS")

#root biomass (most likely perf metric)
NO3biomass <- lmer(Nitrate.C ~ biorat + Level.x + (1|Geno.x),
                data = dat4)
summary(NO3biomass) #plant performance not associated
#weird! There are some genotype-specific relationships with nitrate, and those
#are for genos that have differences in above and belowground biomass allocation
#but these patterns are unrelated to biomass...


# Phosphate ----------------------------------------------------------------------
#Not enough data for phosphate to work

# Lithium ----------------------------------------------------------------------
Li <- lmer(Lithium ~ Soil * Level + scale(Dmass_tot) + (1|block) + (1|Geno),
                data = dat1)

check_model(Li)

summary(Li)
#Lithium is bimodal. Best to run using a GAM

# Sodium -----------------------------------------------------------------------
Sodium <- lmer(log(Sodium) ~ Soil * Level + scale(Dmass_tot)+ Geno +
                 + (1|block),
           data = filter(dat1, cat_batch == 1))

check_model(Sodium)

summary(Sodium)
# For Live plants and microbe pots, Sodium increases as soil moisture increases 
# (L, 2-1 = 0.0857, 0.272, 0.002; N, 2-1 = 1.76, 0.660, 0.009). 
# Microbes lead to an increase in sodium with increasing soil moisture. 
#For C pots, sodium decreases at the highest soil
# moisture level (3-2: -1.24, 0.610, 0.044).
# At this increased soil moisture, sterile plants have
# less sodium than live plants (-0.630, 0.322, 0.053) and microbe pots
# (-1.53, 0.693, 0.029). 
saveRDS(Sodium, file = "04_analyses/02_models/Output/NaSoil.RDS")

NaGeno<- lmer(log(Sodium) ~ Soil + Geno * Level + scale(Dmass_tot) + (1|block),
              data = dat2)
check_model(NaGeno)
summary(NaGeno)

saveRDS(NaGeno, file = "04_analyses/02_models/Output/NaGeno.RDS")
#same pattern as before. Sus.

Nabiomass <- lmer(log(Sodium) ~ biorat + Level.x + (1|Geno.x),
                   data = dat4)
summary(Nabiomass)
#body size does not predict pattern in pore-water sodium variation
# Potassium --------------------------------------------------------------------
Potassium <- lmer(log(Potassium) ~ Soil * Level + scale(Dmass_tot) + Geno +
                    (1|block) ,
               data = filter(dat1, cat_batch == 1))

check_model(Potassium)

summary(Potassium)
# After threshold, both N and L treatments have higher potassium than S
# (L, 2-1 = 0.505, 0.246, 0.043; N, 2-1 = 1.34, 0.528, 0.013). 
#C is not different from anything, but has lower potassium at highest moisture
#(3-2 = -0.859, 0.464, 0.067).
# After threshold, L has higher potassium (0.531, 0.207, 0.012). N also increases
# (1.363, 0.502, 0.008). 
saveRDS(Potassium, file = "04_analyses/02_models/Output/KSoil.RDS")


KGeno<- lmer(log(Potassium) ~ Soil + Geno * Level + scale(Dmass_tot) + (1|block),
            data = dat2)
check_model(KGeno)
summary(KGeno)
#Who knows?

saveRDS(KGeno, file = "04_analyses/02_models/Output/KGeno.RDS")

Kbiomass <- lmer(log(Potassium) ~ massB+ Level.x + (1|Geno.x),
                  data = dat4)
summary(Kbiomass)
#ratio slightly predictive
# Magnesium --------------------------------------------------------------------
Magnesium <- lmer(log(Magnesium) ~ Soil * Level + scale(Dmass_tot) + Geno+
                    (1|block),
                  data = filter(dat1, cat_batch == 1))

check_model(Magnesium)

summary(Magnesium)
#C has lower Mg than N (-1.54, 0.734, 0.0387), this is exacerbated at increased
# soil moisture (2-1 = -3.13, 1.76, 0.078). L has higher Mg than C at increased
# soil moisture too (2-1 = 2.35, 1.385, 0.093). S has lower Mg than L 
# (-0.368, 0.197, 0.065). 

saveRDS(Magnesium, file = "04_analyses/02_models/Output/MgSoil.RDS")

MgGeno<- lmer(log(Magnesium) ~ Soil + Geno *Level + scale(Dmass_tot) + (1|block),
           data = dat2)
check_model(MgGeno)
summary(MgGeno)


saveRDS(MgGeno, file = "04_analyses/02_models/Output/MgGeno.RDS")

Mgbiomass <- lmer(log(Magnesium) ~ massB+ Level.x + (1|Geno.x),
                 data = dat4)
summary(Mgbiomass)
#Mg is higher when root:shoot is higher. Smaller roots may = lower Mg
# Calcium ----------------------------------------------------------------------
Calcium <- lmer(log(Calcium) ~ Soil * Level + scale(Dmass_tot) + Geno +
                  (1|block),
                  data = filter(dat1, cat_batch == 1))

check_model(Calcium)

summary(Calcium)
#N has more Ca than C (1.37, 0.822, 0.099). Ca decreases for C at highest moisture
# (-2.21, 1.066, 0.041). It's also at highest soil moisture that S and L
# (planted treatments) have higher Ca than C, (S, 3-2 = 2.978, 1.13, 0.010;
# L, 3-2 = 2.52, 1.11, 0.026). 
# S increases Ca at highest soil moisture (0.766, 0.417, 0.069). 

saveRDS(Calcium, file = "04_analyses/02_models/Output/CaSoil.RDS")

CaGeno<- lmer(log(Calcium) ~ Soil + Geno * Level + scale(Dmass_tot) + (1|block),
            data = dat2)
check_model(CaGeno)
summary(CaGeno)


saveRDS(CaGeno, file = "04_analyses/02_models/Output/CaGeno.RDS")

Cabiomass <- lmer(log(Calcium) ~ biorat+ Level.x + (1|Geno.x),
                  data = dat4)
summary(Cabiomass)
#no trend with biomass

# Total Carbon -----------------------------------------------------------------

Carbon <- lmer(TC_ppm ~ Soil * Level + scale(Dmass_tot) + TC_batch +
                 (1|Geno) + (1|block),
                  data = dat1)

check_model(Carbon)

summary(Carbon)
#Live soil, planted has more carbon in the soil than sterile planted, however
# at soil moisture level #3, live soil planted has ~= carbon to sterile planted.
# C or N plants aren't in this dataset, so it's unknown whether this is a microbe
# dominated effect. Must get Matt to run the pore-water samples for the C and N
# pots. 

saveRDS(Carbon, file = "04_analyses/02_models/Output/CSoil.RDS")

CGeno<- lmer(Carbon ~ Soil + Geno * Level + scale(Dmass_tot) + (1|block),
            data = dat2)
check_model(CGeno)
summary(CGeno)
# Microbes increase carbon (-150.917, 67.99, 0.031)
# There is genotypic variation.

saveRDS(CGeno, file = "04_analyses/02_models/Output/CGeno.RDS")

Cbiomass <- lmer(Carbon ~ biorat+ Level.x + (1|Geno.x),
                  data = dat4)
summary(Cbiomass)
#body size does not predict pore-water carbon.
# Inorganic Carbon (NPOC) ------------------------------------------------------
dat1$Level <- as.factor(dat1$Level)
contrasts(dat1$Level) <- contr.sum(3)
dat1$batch <- as.factor(dat1$batch)
contrasts(dat1$batch) <- contr.sum(3)

dat2$Level <- as.factor(dat2$Level)
contrasts(dat2$Level) <- contr.sum(3)
dat2$batch <- as.factor(dat2$batch)
contrasts(dat2$batch) <- contr.sum(3)

IC <- lmer(NPOC ~ Soil * Level + scale(Dmass_tot) + batch +
                 (1|Geno) + (1|block),
               data = dat1)

check_model(IC)

summary(IC) #nothing to see here...

ICGeno<- lmer(NPOC ~ Soil + Geno*Level + scale(Dmass_tot) + (1|block),
             data = dat2)
check_model(ICGeno)
summary(ICGeno) #add batch to model as a covariate for actual publication

saveRDS(ICGeno, file = "04_analyses/02_models/Output/ICGeno.RDS")
        
# Performance and PCA geochemistry ---------------------------------------------

#Set1, biomass: above, below, and ratio
mod1 <- lmer(log(massA) ~ (Dim.1 + Dim.2)*Geno + Soil+Level +
               (1|table),
             data = perf2)
             
check_model(mod1)

summary(mod1)
#aboveground biomass not related to geochemistry clusters

mod2 <- lmer(log(massB) ~ (Dim.1 + Dim.2)*Geno + Soil + Level +
               (1|table),
             data = perf2)

check_model(mod2)

summary(mod2)
#belowground biomass not related to geochemistry clusters
mod3 <- lmer(log(biorat) ~ (Dim.1 + Dim.2)*Geno + Soil + Level +
               + (1|table),
             data = perf2)

check_model(mod3)

summary(mod3)#nothing to report here

#Set2, gas exchange: PSN and TPN
mod4 <- lmer(PSN ~ Dim.1 + Dim.2 + Soil*Level +
               (1|Geno) + (1|table),
             data = perf1)

check_model(mod4)

summary(mod4)
#no relationships

mod5 <- lmer(TPN ~ Dim.1 + Dim.2 + Soil*Level +
               (1|Geno) + (1|table),
             data = perf1)

check_model(mod5)

summary(mod5)
#no relationships

#RGR

mod6 <- lmer(RGR_d ~ Dim.1 + Dim.2 + Soil*Level +
               (1|Geno) + (1|table),
             data = perf1)

check_model(mod6)

summary(mod6)
#RGR decreases as Dim2 increases (B = -0.004, SE = 0.002, p = 0.028)
