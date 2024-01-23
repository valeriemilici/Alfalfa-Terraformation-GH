### Initialize workspace--------------------------------------------------------
rm(list =ls())
library(lme4)
library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

ICmod<- readRDS("04_analyses/02_models/Output/ICSoil.RDS")
load("04_Analyses/03_Bootstrap/output/IC.boot")


modpredIC <- function(mod, boot){
  pred <- expand.grid(Soil = c("C", "N", "S", "L"),
                      Level = c("1", "2", "3"),
                      Dmass_tot = 106,
                      batch = c("A", "B", "C"))
  pred$preds <- predict(mod, newdata = pred, re.form = ~0)
  pred <- data.frame(pred, confint(boot))
  names(pred)[6:7] <- c("lwr", "upr")
  return(pred)
}

predIC <- modpredIC(ICmod, IC.boot)


# Create Plot

soil.labs <- c("No Microbe | No Plant", "Yes Microbes | No Plant",
               "No Microbes | Yes Plant", "Yes Microbes | Yes Plant")
names(soil.labs) <- c("C", "N", "S", "L")


ICplot <- ggplot() +
  geom_pointrange(data = filter(predIC, batch == "B"),
                  mapping = aes(Level, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  facet_wrap(~Soil, labeller = labeller(Soil = soil.labs)) + 
  xlab("Soil Moisture Treatment") +
  ylab("Inorganic Carbon (ug/ml)") +
  theme_classic(20) 

IClot
ggsave(ICplot, filename = "figures/IC.png")
