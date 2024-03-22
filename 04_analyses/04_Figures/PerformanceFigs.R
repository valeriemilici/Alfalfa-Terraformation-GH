### Make some figures!

### Initialize workspace--------------------------------------------------------
rm(list =ls())
library(lme4)
library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(patchwork) #joins plots nicely
library(ggeffects) #plot directly from model without bootstrap

### Load Bootstrap Output ------------------------------------------------------

load("04_analyses/03_Bootstrap/output/TotBM.boot")
load("04_analyses/03_Bootstrap/output/RMF.boot")
#load("04_analyses/03_Bootstrap/output/Mort.boot")

### Load Models ----------------------------------------------------------------

TotBM <- readRDS("04_analyses/02_models/Output/Biomass.RDS")
RMF <- readRDS("04_analyses/02_models/Output/RMF.RDS")
RGR <- readRDS("04_analyses/02_models/Output/RGR.RDS")
Mort <- readRDS("04_analyses/02_models/Output/Mortality.RDS")

### General Biomass Trends -------------------------------------------------

predBiomass <- expand.grid(Soil = c("L", "S"),
                       Level = c("1", "2", "3"))

predBiomass$preds <- predict(TotBM, newdata = predBiomass, re.form = ~0)

predBiomass <- data.frame(predBiomass, confint(TotBM.boot))

names(predBiomass)[4:5] <- c("lwr", "upr")

# Create Plot

Biomassplot <- ggplot() +
  geom_pointrange(data = predBiomass,
                  mapping = aes(Level, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr),
                                group = Soil, col = Soil),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Percent of Water Holding Capacity") +
  ylab("Total Dry Biomass (g)") +
  scale_x_discrete(labels = c("1" = "30-50", "2" = "55-75", "3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(15)

Biomassplot

ggsave(Biomassplot, filename = "figures/Biomass.png",
       height = 6, width =8, units = "in")

### Root Mass Fraction Trends -------------------------------------------------

predRMF <- expand.grid(totBM = 0.089,
                       Soil = c("L", "S"),
                           Level = c("1", "2", "3")
                       )

predRMF$preds <- predict(RMF, newdata = predRMF, re.form = ~0)

predRMF <- data.frame(predRMF, confint(RMF.boot))

names(predRMF)[5:6] <- c("lwr", "upr")

# Create Plot

RMFplot <- ggplot() +
  geom_pointrange(data = predRMF,
                  mapping = aes(Level, preds,
                                ymin = lwr, ymax = upr,
                                group = Soil, col = Soil),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Percent of Water Holding Capacity") +
  ylab("Root Mass Fraction") +
  scale_x_discrete(labels = c("1" = "30-50", "2" = "55-75", "3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(15)

RMFplot

ggsave(RMFplot, filename = "figures/RMF.png",
       height = 6, width =8, units = "in")

### Relative Growth Rate Trends -------------------------------------------------

#try this out with ggeffects

rgrResult <- ggemmeans(RGR, c("Level", "Soil"), type = "fixed")

# Create Plot

RGRplot <- ggplot() +
  geom_pointrange(data = rgrResult,
                  mapping = aes(x, predicted,
                                ymin = conf.low, ymax = conf.high,
                                group = group, col = group),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Percent of Water Holding Capacity") +
  ylab("Relative Growth Rate (mm/mm/day)") +
  scale_x_discrete(labels = c("1" = "30-50", "2" = "55-75", "3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(15)

RGRplot

ggsave(RGRplot, filename = "figures/RGR.png",
       height = 6, width =8, units = "in")

### Mortality Trends -----------------------------------------------------------

MortPreds <- ggpredict(Mort, c("Level", "Soil"))

Mortplot <- ggplot() +
  geom_pointrange(data = MortPreds,
                  mapping = aes(x, predicted,
                                ymin = conf.low, ymax = conf.high,
                                group = group, col = group),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Percent of Water Holding Capacity") +
  ylab("Probability of Mortality") +
  scale_x_discrete(labels = c("1" = "30-50", "2" = "55-75", "3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(15)

Mortplot

ggsave(Mortplot, filename = "figures/Mort.png",
       height = 6, width =8, units = "in")

### Genotype-Specific Biomass Ratio --------------------------------------------

# Predict new data & add CIs. Creates data for plotting.
predBR <- expand.grid(Soil = "L",
                      totBM = 0.067,
                          Level = c("1","2","3"),
                          Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))

predBR$preds <- predict(BR, newdata = predBR, re.form = ~0)

predBR<- data.frame(predBR, confint(BR.boot))

names(predBR)[6:7] <- c("lwr", "upr")

# Create Plot

BRplot1 <- ggplot() +
  geom_pointrange(data = filter(predBR, Level == "1"),
                  mapping = aes(Geno, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr)),
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = exp(fixef(BR)[1] + fixef(BR)[2]*log(0.067)), col = "red") +
  xlab("") +
  ylab("Root:Shoot Mass (g)") +
  ylim(0,1) +
  theme_classic(12) 

BRplot1

BRplot2 <- ggplot() +
  geom_pointrange(data = filter(predBR, Level == "2"),
                  mapping = aes(Geno, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr)),
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = exp(fixef(BR)[1] + fixef(BR)[4]+
                                fixef(BR)[2]*log(0.067)), col = "red") +
  xlab("Genotypes") +
  ylab("Root:Shoot Mass (g)") +
  ylim(0,1) +
  theme_classic(12) 

BRplot2

BRplot3 <- ggplot() +
  geom_pointrange(data = filter(predBR, Level == "3"),
                  mapping = aes(Geno, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr)),
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = exp(fixef(BR)[1] +fixef(BR)[5]+
                                fixef(BR)[2]*log(0.067)), col = "red") +
  xlab("") +
  ylab("") +
  ylim(0,1) +
  theme_classic(12)

BRplot3

BR.Moisture <- BRplot1 +
  BRplot2 + 
  BRplot3 +
  plot_annotation(caption = "Alfalfa Genotype",
                  theme = theme(plot.caption = element_text(size = 22,
                                                            hjust = 0.5)))

BR.Moisture

ggsave(BRplot2, filename = "figures/BR.Moisture2.png")
### Genotype Specific Total Biomass --------------------------------------------

# Predict new data & add CIs. Creates data for plotting.
predBT <- expand.grid(Soil = c("L", "S"),
                      Level = c("2"),
                      Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))

predBT$preds <- predict(BT, newdata = predBT, re.form = ~0)

predBT<- data.frame(predBT, confint(BT.boot))

names(predBT)[5:6] <- c("lwr", "upr")

# Create Plot

BTplot <- ggplot() +
  geom_pointrange(data = filter(predBT, Soil == "L"),
                  mapping = aes(Geno, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr)),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = exp(fixef(BT)[1]), col = "red") +
  xlab("Alfalfa Genotypes") +
  ylab("Total Biomass (g)") +
  theme_classic() +
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text  = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

BTplot
ggsave(BTplot, filename = "figures/BT.jpeg",
       height = 8, width =10, units = "in")
### Genotype Specific Photosynthesis -------------------------------------------

predPSN <- expand.grid(Soil = c("L", "S"),
                      Level = c("2"),
                      Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))

predPSN$preds <- predict(PSN, newdata = predPSN, re.form = ~0)

predPSN<- data.frame(predPSN, confint(PSN.boot))

names(predPSN)[5:6] <- c("lwr", "upr")

# Create Plot

PSNplot <- ggplot() +
  geom_pointrange(data = filter(predPSN, Soil == "L"),
                  mapping = aes(Geno, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = fixef(PSN)[1], col = "red") +
  xlab("Alfalfa Genotypes") +
  ylab("Photosynthesis \n(Scott, what are the units?)") +
  theme_classic() +
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text  = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

PSNplot

ggsave(PSNplot, filename = "figures/PSN.jpeg",
       height = 8, width =10, units = "in")
### Genotype Specific Transpiration -------------------------------------------

predTPN <- expand.grid(Soil = c("L", "S"),
                       Level = c("2"),
                       Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))

predTPN$preds <- predict(TPN, newdata = predTPN, re.form = ~0)

predTPN<- data.frame(predTPN, confint(TPN.boot))

names(predTPN)[5:6] <- c("lwr", "upr")

# Create Plot

TPNplot <- ggplot() +
  geom_pointrange(data = filter(predTPN, Soil == "L"),
                  mapping = aes(Geno, preds,
                                ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = fixef(TPN)[1], col = "red") +
  xlab("Genotypes") +
  ylab("Transpiration \n (mmol/m^2/s)") +
  theme_classic(12) 

TPNplot

ggsave(TPNplot, filename = "figures/TPN.png")

### Genotype Specific Survival -------------------------------------------------

predSurv <- expand.grid(Soil = c("L", "S"),
                       Level = c("2"),
                       Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))

predSurv$preds <- predict(Surv, newdata = predSurv, re.form = ~0)

predSurv <- data.frame(predSurv, confint(Surv.boot))

names(predSurv)[5:6] <- c("lwr", "upr")

# Create Plot

Survplot <- ggplot() +
  geom_pointrange(data = filter(predSurv, Soil == "L"),
                  mapping = aes(Geno, plogis(preds),
                                ymin = plogis(lwr), ymax = plogis(upr)),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = plogis(fixef(Surv)[1]), col = "red") +
  xlab("Alfalfa Genotypes") +
  ylab("Prob. of Survival") +
  theme_classic() +
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text  = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

Survplot
ggsave(Survplot, filename = "figures/Surv.jpeg",
       height = 8, width =10, units = "in")

