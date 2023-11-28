## Models specifically on root traits and performance

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(MASS) #for contr.sdif
library(performance) #model diagnostics

dat <- read.csv("data/performance.csv")
alfdat <- read.csv("data/GenoSelections/AlfPISelections.csv") #alfalfa data

### clean the data--------------------------------------------------------------
#Remove NAs and simplify data to one observation per plant
dat1 <- dat %>% filter(massA > 0) %>%
  group_by(pot) %>%
  arrange(CensusNo) %>%
  #filters out so you have the last census for each pot
  filter(row_number() == n()) %>%
  mutate(totBM = massA + massB,
         #center the responses to make them more "normal" w/o logging
         totBM.c = totBM - mean(totBM),
         massA.c = massA - mean(massA),
         massB.c = massB - mean(massB))

#Prep alf dat so that it has geno column/codes, select the relevent columns
Geno <- c("D16", "A20", "G15", "212", "SWE" ,"R24", "TAS", "YON", "K19",
          "MOC", "INA", "CRI", "639", "105", "G74", "515", "VIR")
alfdat2 <- cbind(Geno, alfdat)

alfdat3<- alfdat2 %>% select(Geno, Days.Plant.to.Maturity,
                             Determinate.Taproot.Percentage,
                             Determinate.Taproot.Position, Fibrous.Root.Mass,
                             Secondary.Root.Number)

#bring them together!
dat2 <- left_join(dat1, alfdat3, by = "Geno") 

dat3 <- dat2 %>% mutate(DTP = ifelse(Determinate.Taproot.Percentage <=40,
                                     "A", "B"))

#manipulations for the survival models ---------
dat$HtPrevCensus[dat$CensusNo == 1] <- 0 #set Census 1 prev ht to 0 not NA

dat4 <- dat %>%
  filter(HtPrevCensus >= 0) %>% 
  group_by(pot) %>% 
  arrange(CensusNo) %>%
  filter(row_number() == 1 | row_number() == n())

dat5 <- left_join(dat4, alfdat3, by = "Geno")
dat6 <- dat5 %>% mutate(DTP = ifelse(Determinate.Taproot.Percentage <=40,
                                      "A", "B"))

#why is this here? Does "probs" mean problems? Or Probability?
#probs <- dat1 %>% group_by(pot) %>%
 # summarise(n = n()) %>%
  #filter(n == 1)

### The models -----------------------------------------------------------------

##set contrasts-----------------------------------------------------------------

dat3$Level <- as.factor(dat3$Level)
contrasts(dat3$Level) <- contr.sdif(3)

## [1]Biomass Mods-------------------------------

# Total Biomass ----------------
mod1 <- lmer(log(totBM) ~ (DTP +
               Fibrous.Root.Mass) * Soil + Level +
               (1|Geno) + (1|table),
                data = dat3)

check_model(mod1)
summary(mod1)
#In sterile soil, plants with a higher fibrous root mass had a larger overall
#biomass. A sign that something was attacking the roots in the live soil?

# Aboveground Biomass --------
mod2 <- lmer(log(massA) ~ (DTP +
                             Fibrous.Root.Mass) * Soil + Level +
               (1|Geno) + (1|table),
             data = dat3)

check_model(mod2)
summary(mod2)
#It seems like the difference in biomass when comparing sterile to live plants
#is owed to sterile soil plants with fibrous genotypes having a larger root
#system.

# Belowground Biomass ------
mod3a <- lmer(log(massB) ~ log(totBM) + (DTP +
                             Fibrous.Root.Mass) + Soil + Level +
               (1|Geno) + (1|table),
             data = dat3,
             REML = F)

mod3b <- lmer(log(massB) ~ (DTP +
                                           Fibrous.Root.Mass) * Soil + Level +
                (1|Geno) + (1|table),
              data = dat3,
              REML = F)
anova(mod3a, mod3b)
check_model(mod3a)
summary(mod3a)
#revision! The results completely change when I balance out the model with 
#total biomass on the RHS. This model is much better (via anova) than the OG
#model without biomass. And now it shows that all patterns are owed to bigger
#overall plants having bigger roots, and there is no relationship with root
#traits to explain this pattern. Now, the bigger plants were the sterile plants
#and that is not coming out as anything anymore, but those two 

# Root to shoot ratio -------

mod3.1<- lmer(log(massA/massB) ~  log(totBM) +
                (DTP + Fibrous.Root.Mass) * Soil + Level +
                (1|Geno) + (1|table),
               data = filter(dat3, rootsample == "Y"),
              REML = F)

summary(mod3.1)

##[2] Gas Exchange Mods -----------------------

# PSN --------------
mod4 <- lmer(PSN ~ (DTP +
                      Fibrous.Root.Mass) * Soil + Level +
               (1|Geno) + (1|table),
             data = filter(dat3))

check_model(mod4)
# this model looks better when all of the data are included, instead of filtered
# such that PSN < 0.07. However, this still isn't perfect. The model seems to 
# underpredict the high values, and hte PPC is only a decent match for the
# observed data. 

summary(mod4)
# Nothing much happening here, which is ok with me. We just find that PSN 
# increases between level 1 and 2, which is consistent with all of our other
# findings that the biggest benefit for plants was between level 1 and 2, and 
# that performance in 2&3 is pretty similar. 

# TPN -------------------
mod5 <- lmer(TPN ~ (DTP +
                      Fibrous.Root.Mass) * Soil + Level +
               (1|Geno) + (1|table),
             data = filter(dat3))

check_model(mod5)
# The model does a bad job of predicting the high values, otherwise it looks
# pretty good.

summary(mod5)
# nothing really going on here with transpiration rates. 

## Overall it seems like the root traits that we EXPECT the plants to have don't
## have any relationship with their gas exchange realities during the experiment.

##[3] Survival ------------------------------

mod6 <- glmer(as.numeric(status) ~  Fibrous.Root.Mass  + Soil + Level +
                      (1|Geno) + (1|table),
                      data = dat6,
                      family = binomial(link = "logit"))

check_model(mod6) #not sure if the performance package is the best
#diagnostic package


summary(mod6)
#In general nothing other than soil moisture and maybe sterile soil affects
#survival. Root traits aren't associated with survival. 
