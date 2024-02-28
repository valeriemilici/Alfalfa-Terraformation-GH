### Preliminary Analysis of Soil Sample Geochemistry on Priority Samples Only

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
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

dat <- read.csv("data/SolidSamples_C_N_PriorityOnly.csv") #geochem data
location <- read.csv("data/GCR_2A_CensusData.csv") #sample metadata

### Prepare data for models ----------------------------------------------------
location1 <- location %>% 
  filter(CensusNo == "6") %>%
  dplyr::select(c(pot, table, rack, position)) #gets only what we need

dat1 <- dat %>%
  rename(pot = Sample.ID) %>%
  group_by(pot, Acid.Treatment) %>%
  summarize(Carbon = mean(TC_ug.mg),
            Nitrogen = mean(TN_ug.mg),
            batch = mean(batch)) #enables join

dat2 <- left_join(dat1, location1, by = "pot") %>% #combine data
  #create new columns for model
  mutate(Soil = str_extract(pot, "\\w"),
         Geno = str_extract(pot, "\\w{2,}"),
         Level = substr(pot, 3,3),
         Location = str_c(table, rack, sep = ".")) 

#separate data based on two different methods
datY <- dat2 %>% filter(Acid.Treatment == "YES")
datN <- dat2 %>% filter(Acid.Treatment == "NO")

### Set Contrasts --------------------------------------------------------------

datY$Level <- as.factor(datY$Level)
contrasts(datY$Level) <- contr.sdif(3) 
#so models can be interpreted as step-wise increases in water availability

datN$Level <- as.factor(datN$Level)
contrasts(datN$Level) <- contr.sdif(3)
#so models can be interpreted as step-wise increases in water availability

### Relevel "Soil" treatment to understand patterns ----------------------------
datY$Soil <- relevel(as.factor(datY$Soil),
                     ref = "L") 

datN$Soil <- relevel(as.factor(datN$Soil),
                     ref = "N") #Switch between "L", "C", and "N"

### Model 1: Patterns with total Carbon ----------------------------------------

mod1 <- lmer(Carbon ~ Soil + Level + (1|Geno) + (1|batch) + (1|Location), 
             data = datN)

check_model(mod1)
anova(mod1)

#confirmed through model ANOVA and AIC comparison of two models that Soil and 
#level should not interact in the model. 

summary(mod1)

### Create Figure of Model 1 ---------------------------------------------------
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

# The functions -----------------
geochem.fun <- function(.) {
  preddat <- expand.grid(Soil = c("C", "N", "S", "L"),
                         Level =c("1", "2", "3"))
  predict(., newdata = preddat, re.form = ~0)
}

modpred <- function(mod, boot){
  pred <- expand.grid(Soil = c("C", "N", "S", "L"),
                      Level = c("1", "2", "3"))
  pred$preds <- predict(mod, newdata = pred, re.form = ~0)
  pred <- data.frame(pred, confint(boot))
  names(pred)[4:5] <- c("lwr", "upr")
  return(pred)
}

# Create the model predictions to plot ----
Cboot <- bootMer(mod1, nsim = 1000, FUN = geochem.fun,
                      parallel="snow", ncpus = detectCores(), 
                      cl=cl) #bootstrap the model (to get CIs)

stopCluster(cl = cl) 

predC <- modpred(mod1, Cboot)

# Build the plot ----
soil.labs <- c("No Microbe | No Plant", "Yes Microbes | No Plant",
               "No Microbes | Yes Plant", "Yes Microbes | Yes Plant")
names(soil.labs) <- c("C", "N", "S", "L")
               
plot1 <- ggplot() +
  geom_pointrange(data = predC,
                  mapping = aes(Level, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  facet_wrap(~Soil, labeller = labeller(Soil = soil.labs)) + 
  xlab("Soil Moisture Treatment") +
  ylab("Total Carbon (ug/mg)") +
  theme_classic(20)

plot1

ggsave(plot = plot1, filename = "figures/SolidPhaseCarbon.png")
