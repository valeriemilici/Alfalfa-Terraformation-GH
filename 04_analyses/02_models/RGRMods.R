### This is an analysis on the 2A alfalfa data that asks...
### is there a difference in  RGR between seedlings growing in live vs. sterile
### & by moisture/genos?

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

perfdat <- dat %>% filter(!Geno %in% outs & Status == "A" &
                            CensusNo != 1) 
#RGR only possible with live plants & at least 2 censuses of data

# 5 rows appear to be outliers
dat2 <- dat[-c(680, 1721, 1260, 679, 244), ]
perfdat2 <- perfdat[-c(207,63,510,26,208),]
### setting contrasts ----------------------------------------------------------

perfdat2$Geno <- as.factor(perfdat2$Geno)
perfdat2$Geno <- factor(perfdat2$Geno,
                       levels = c("CRI", "INA", "K19", "MOC",
                                  "TAS", "VIR", "YON", "G15"))
perfdat$Level <- as.factor(perfdat$Level)
contrasts(perfdat2$Geno) <- contr.sum(8) #intercept represents ave across genos
contrasts(perfdat$Level) <- contr.helmert(3) #intercept is ave across soil water

dat2$Level <- as.factor(dat2$Level)
contrasts(dat2$Level) <- contr.helmert(3) #this reduces collinearity issue
### models ---------------------------------------------------------------------

## I am worried about having enough data to really ask these questions. The 
## models will increase in complexity as I go forward, so I can learn where to 
## stop.

### model 1: Soil*Level across the entire experiment ---------------------------
dat2$Soil <- factor(dat2$Soil, levels = c("S", "L"))

RGRsoil <- lmer(RGR_d ~ Soil*Level +
                  HtPrevCensus +
                  (1|CensusNo) + (1|pot) + (1|Geno) + (1|table),
                data = filter(dat2, CensusNo != 1 & Status == "A"))

## model 1 diagnostics -----------------------------------
#[1] Summary of model diagnostics
check_model(RGRsoil)
# removal of 5 points has improved q-q plot
#[2] Influential outliers?
ggplot(data.frame(lev=hatvalues(RGRsoil),pearson=residuals(RGRsoil,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw() 

levId <- which(hatvalues(RGRsoil) >= .03) #rows w/potentially high influence
#the model is pointing out L-3-YON-06 as being weird, and it goes from very
#high RGR to very low RGR
dat3 <- filter(dat2, pot != "L-3-YON-06")
#WTH this isn't right! I keep getting the same number of values that have the 
#high levIDs! This is upsetting and I think I need to wait a moment

# model 1 output -----------------------------------
summary(RGRsoil) 

#across the duration of the entire experiment, plants in sterile soil grew
#faster than plants in live soil (B = 6.4e-3, std =  1.19e-3, P < 0.001). 
#all plants grew faster as soil moisture increased and this benefit was 
# regardless of whether the plant was in live or sterile soil. 

### model 2: Soil*Level at the beginning of the experiment----------------------

RGRsoil2 <- lmer(RGR_d ~ Soil*Level +
                  HtPrevCensus +
                  (1|Geno) + (1|table),
                data = filter(dat2, CensusNo == 2 & Status == "A"))

## model 2 diagnostics -----------------------------------
#[1] Summary of model diagnostics
check_model(RGRsoil2) #looks good!
#[2] Influential outliers?
ggplot(data.frame(lev=hatvalues(RGRsoil2),
                  pearson=residuals(RGRsoil2,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw() 
#this looks like a good spread of points!

# model 2 output -----------------------------------
summary(RGRsoil2) 

#when we only consider the initial censuses, early in the experiment, we see the
#same trend as above, except that there is a marginal interaction between 
#sterile soil and soil moisture, which tells us that early on the RGR of
# plants in sterile soil were not advantaged by an increase in soil moisture, 
# and this advantage was only experienced by plants in Live soil. 

### model 3: Genotype-specific responses ---------------------------------------
RGRgeno <- lmer(RGR_d ~ Soil + Geno*Level +
                   HtPrevCensus +
                   (1|pot) + (1|table) + (1|CensusNo),
                 data = perfdat2)

#ggplot(data.frame(lev=hatvalues(RGRgeno),
#                  pearson=residuals(RGRgeno,type="pearson")),
#       aes(x=lev,y=pearson)) +
#  geom_point() +
#  theme_bw() 

#not necessary to see, but if you're curious about my cut-off and want to play 
# with others, go ahead.

levId <- which(hatvalues(RGRgeno) >= .06) #15 data points

RGRgeno1 <- lmer(RGR_d ~ Soil + Geno*Level +
                   HtPrevCensus +
                   (1|pot) + (1|table) + (1|CensusNo),
                 data = perfdat2[-c(as.vector(levId)), ]) 
#model without those high hat values

## model 3 diagnostics -----------------------------------
#[1] Summary of model diagnostics
check_model(RGRgeno1) #improved from before

# model 3 output -----------------------------------
summary(RGRgeno1)  

# CRI genotype has marginally higher RGR in LEO soil than
# other genotypes (p = 0.08). No sig interactions but CRI has a curious
# almost-marginal where it grows faster as moisture increases
# (p = 0.13). Considering the sample size, that is compelling. 
# It's interesting that this fast overall RGR doesn't make it a better survivor

### model 4: Genotype-specific early on ----------------------------------------
RGRgeno2 <- lmer(RGR_d ~ Soil + Geno*Level +
                  HtPrevCensus +
                  (1|table),
                data = filter(perfdat2, CensusNo == 2))
## model 4 diagnostics ----------------
check_model(RGRgeno2) #decent

summary(RGRgeno2)
# early on, YON had marginally above-average RGR (p = 0.09).
# INA grew slowly at intermediate moisture (p = 0.03) as did VIR (p = 0.09).
# K19 grew quickly at the highest soil moisture (p = 0.09). 
