### Make some figures!

### Initialize workspace--------------------------------------------------------
rm(list =ls())
library(lme4)
library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

### Load Bootstrap Output ------------------------------------------------------

load("04_analyses/03_Bootstrap/output/BR.boot")
load("04_analyses/03_Bootstrap/output/BT.boot")
load("04_analyses/03_Bootstrap/output/BSR.boot")
load("04_analyses/03_Bootstrap/output/BST.boot")
load("04_analyses/03_Bootstrap/output/PSN.boot")
load("04_analyses/03_Bootstrap/output/TPN.boot")
load("04_analyses/03_Bootstrap/output/Surv.boot")

### Load Models ----------------------------------------------------------------

BR <- readRDS("04_analyses/02_models/output/BMGenoR.RDS")
BT <- readRDS("04_analyses/02_models/output/BMGenoT.RDS")
BST <- readRDS("04_analyses/02_models/output/BMSoilT.RDS")
BSR <- readRDS("04_analyses/02_models/output/BMSoilR.RDS")
PSN <- readRDS("04_analyses/02_models/output/PSNGeno.RDS")
TPN <- readRDS("04_analyses/02_models/output/TPNGeno.RDS")
Surv <- readRDS("04_analyses/02_models/output/SurvGeno1.RDS")


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

### Biomass Geno Blind -------------------------------------------------

predBST <- expand.grid(Soil = c("L", "S"),
                        Level = c("1", "2", "3"))

predBST$preds <- predict(BST, newdata = predBST, re.form = ~0)

predBST <- data.frame(predBST, confint(BST.boot))

names(predBST)[4:5] <- c("lwr", "upr")

# Create Plot

BSTplot <- ggplot() +
  geom_pointrange(data = predBST,
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
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

BSTplot
ggsave(BSTplot, filename = "figures/BST.jpeg",
       height = 6, width =8, units = "in")

