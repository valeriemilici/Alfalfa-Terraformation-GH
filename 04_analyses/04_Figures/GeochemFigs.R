### Make some figures!

### Initialize workspace--------------------------------------------------------
rm(list =ls())
library(lme4)
library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(patchwork) #joins plots nicely

### Load Bootstrap Output ------------------------------------------------------
load("04_analyses/03_Bootstrap/output/EC.boot")
load("04_analyses/03_Bootstrap/output/pH.boot")
load("04_analyses/03_Bootstrap/output/C.boot")
load("04_analyses/03_Bootstrap/output/Ca.boot")
load("04_analyses/03_Bootstrap/output/K.boot")
load("04_analyses/03_Bootstrap/output/Mg.boot")
load("04_analyses/03_Bootstrap/output/Na.boot")
load("04_analyses/03_Bootstrap/output/NO3.boot")

load("04_Analyses/03_Bootstrap/output/ECGeno.boot")
load("04_Analyses/03_Bootstrap/output/CGeno.boot")
load("04_Analyses/03_Bootstrap/output/CaGeno.boot")
load("04_Analyses/03_Bootstrap/output/MgGeno.boot")
load("04_Analyses/03_Bootstrap/output/KGeno.boot")
load("04_Analyses/03_Bootstrap/output/NaGeno.boot")
load("04_Analyses/03_Bootstrap/output/NO3Geno.boot")
load("04_Analyses/03_Bootstrap/output/ICGeno.boot")
### Load Models ----------------------------------------------------------------
EC <- readRDS("04_analyses/02_models/Output/ECSoil.RDS")
pH <- readRDS("04_analyses/02_models/Output/pHSoil.RDS")
Carbon <- readRDS("04_analyses/02_models/Output/CSoil.RDS")
Calcium <- readRDS("04_analyses/02_models/Output/CaSoil.RDS")
Potassium <- readRDS("04_analyses/02_models/Output/KSoil.RDS")
Magnesium <- readRDS("04_analyses/02_models/Output/MgSoil.RDS")
Sodium <- readRDS("04_analyses/02_models/Output/NaSoil.RDS")
Nitrate <- readRDS("04_analyses/02_models/Output/NO3Soil.RDS")

ECGeno <- readRDS("04_analyses/02_models/Output/ECGeno.RDS")
CGeno <- readRDS("04_analyses/02_models/Output/CGeno.RDS")
CaGeno <- readRDS("04_analyses/02_models/Output/CaGeno.RDS")
KGeno <- readRDS("04_analyses/02_models/Output/KGeno.RDS")
MgGeno <- readRDS("04_analyses/02_models/Output/MgGeno.RDS")
NaGeno <- readRDS("04_analyses/02_models/Output/NaGeno.RDS")
NO3Geno <- readRDS("04_analyses/02_models/Output/NO3Geno.RDS")
ICGeno <- readRDS("04_analyses/02_models/Output/ICGeno.RDS")
### Prediction Matrix Function -------------------------------------------------

modpred <- function(mod, boot){
  pred <- expand.grid(Soil = c("C", "N", "S", "L"),
                      Level = c("1", "2", "3"),
                      Dmass_tot = 106)
  pred$preds <- predict(mod, newdata = pred, re.form = ~0)
  pred <- data.frame(pred, confint(boot))
  names(pred)[5:6] <- c("lwr", "upr")
  return(pred)
}

modpredG <- function(mod,boot){
  pred <- expand.grid(Soil = c("S","L"),
                      Level = c("1","2", "3"),
                      Dmass_tot = 106,
                      Geno = c("CRI", "K19", "MOC",
                              "TAS", "VIR", "G15", "YON"))
  pred$preds <- predict(mod, newdata = pred, re.form = ~0)
  pred <- data.frame(pred, confint(boot))
  names(pred)[6:7] <- c("lwr", "upr")
  pred <- filter(pred, Soil == "L" & Level == "2")
 
  pred <- pred[order(pred$preds),]
  newlevel <- pred$Geno
  pred$Geno <- factor(pred$Geno,
                      levels = newlevel)
  return(pred)
}

modpredC <- function(mod, boot){
  pred <- expand.grid(Soil = c("C", "N", "S", "L"),
                      Level = c("1", "2", "3"),
                      Dmass_tot = 106,
                      TC_batch = "1")
  pred$preds <- predict(mod, newdata = pred, re.form = ~0)
  pred <- data.frame(pred, confint(boot))
  names(pred)[6:7] <- c("lwr", "upr")
  return(pred)
}

modpredCat <- function(mod, boot){
  pred <- expand.grid(Soil = c("S", "L"),
                      Level = c("1", "2", "3"),
                      Dmass_tot = 106,
                      Geno = c("MOC", "VIR", "G15", "105"))
  pred$preds <- predict(mod, newdata = pred, re.form = ~0)
  pred <- data.frame(pred, confint(boot))
  names(pred)[6:7] <- c("lwr", "upr")
  return(pred)
}

  
### Create the Figures ---------------------------------------------------------
## EC Soil ---------------------------------------------------------------------

predEC <- modpred(EC, EC.boot)

# Create Plot

ECplot <- ggplot() +
  geom_pointrange(data = filter(predEC, Level == "1"),
                  mapping = aes(Soil, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr)),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Soil Treatment") +
  ylab("Electrical Conductivity (EC)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

ECplot
ggsave(ECplot, filename = "figures/EC.jpeg",
       height = 6, width =8, units = "in")

## EC Geno ---------------------------------------------------------------------
ECpredG <- modpredG(ECGeno, ECGeno.boot)

ins <- c("VIR", "YON", "TAS", "CRI", "MOC")

ECpredG2 <- ECpredG %>% filter(Geno %in% ins) 

ECGplot <- ggplot() +
  geom_pointrange(data = ECpredG2,
                  mapping = aes(Geno, exp(preds)/1000,
                                ymin = exp(lwr)/1000, ymax = exp(upr)/1000),
                  size = 0.9,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = exp(fixef(ECGeno)[1])/1000, col = "red") +
  xlab("Alfalfa Cultivar") +
  ylab("Electrical Conductivity (mS/m)") +
  theme_classic() +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text  = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10))

ECGplot
ggsave(ECGplot, filename = "figures/ECGeno.jpeg",
       height = 2, width =2.67, units = "in")

## pH Soil ---------------------------------------------------------------------

predpH <- modpred(pH, pH.boot)

# Create Plot

pHplot <- ggplot() +
  geom_pointrange(data = filter(predpH, Level == "2"),
                  mapping = aes(Soil, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Soil Treatment") +
  ylab("pH") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

pHplot
ggsave(pHplot, filename = "figures/pH.jpeg",
       height = 6, width =8, units = "in")
## Carbon Soil -----------------------------------------------------------------


predC <- modpredC(Carbon, C.boot)


# Create Plot

soil.labs <- c("No Microbe | No Plant", "Yes Microbes | No Plant",
               "No Microbes | Yes Plant", "Yes Microbes | Yes Plant")
names(soil.labs) <- c("C", "N", "S", "L")


Carbonplot <- ggplot() +
  geom_pointrange(data = predC,
                  mapping = aes(Level, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  facet_wrap(~Soil, labeller = labeller(Soil = soil.labs)) + 
  xlab("Soil Moisture Treatment") +
  ylab("Total Carbon (ppm)") +
  theme_classic(20) 

Carbonplot
ggsave(Carbonplot, filename = "figures/Carbon.png")

## Carbon Geno -----------------------------------------------------------------
CpredG <- modpredG(CGeno, CGeno.boot)

ins <- c("VIR", "YON", "TAS", "CRI", "MOC")

CpredG2 <- CpredG %>% filter(Geno %in% ins) 

CGplot <- ggplot() +
  geom_pointrange(data = CpredG2,
                  mapping = aes(Geno, preds,
                                ymin = lwr, ymax = upr),
                  size = 0.9,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = fixef(CGeno)[1], col = "red") +
  xlab("Genotype") +
  ylab("Total Carbon (ppm)") +
  theme_classic() +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text  = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10))

CGplot
ggsave(CGplot, filename = "figures/CarbonGeno.jpeg",
       height = 3, width =4, units = "in")
## Inorganic Carbon Geno -------------------------------------------------------
ICpredG <- modpredG(ICGeno, ICGeno.boot)

ins <- c("VIR", "YON", "TAS", "CRI", "MOC")

ICpredG2 <- ICpredG %>% filter(Geno %in% ins) 

ICGplot <- ggplot() +
  geom_pointrange(data = ICpredG2,
                  mapping = aes(Geno, preds,
                                ymin = lwr, ymax = upr),
                  size = 0.9,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = fixef(ICGeno)[1], col = "red") +
  xlab("Alfalfa Cultivar") +
  ylab("Inorganic Carbon (ppm)") +
  theme_classic() +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text  = element_text(size = 10)) +
  theme(legend.text = element_text(size = 10))

ICGplot
ggsave(ICGplot, filename = "figures/InorganicCarbonGeno.jpeg",
       height = 2, width =2.67, units = "in")
## Calcium Soil ----------------------------------------------------------------
predCa <- modpredCat(Calcium, Ca.boot)

# Create Plot

Caplot <- ggplot() +
  geom_pointrange(data = filter(predCa, Level == "3"),
                  mapping = aes(Soil, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr),
                                group = Geno, col = Geno),
                  size = 1.2,
                  position = position_dodge(width = 0.5)) +
  xlab("Soil Treatment") +
  ylab("Calcium (mg/L)") +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

Caplot
ggsave(Caplot, filename = "figures/Ca.jpeg",
       height = 6, width =8, units = "in")

## Calcium Geno ----------------------------------------------------------------
CapredG <- modpredG(CaGeno, CaGeno.boot)

CaGplot <- ggplot() +
  geom_pointrange(data = CapredG,
                  mapping = aes(Geno, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = fixef(CaGeno)[1], col = "red") +
  xlab("Genotype") +
  ylab("log Calcium (mg/L)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

CaGplot
ggsave(CaGplot, filename = "figures/CalciumGeno.jpeg",
       height = 6, width =8, units = "in")
## Potassium Soil --------------------------------------------------------------
predK <- modpredCat(Potassium, K.boot)

# Create Plot

Kplot <- ggplot() +
  geom_pointrange(data = filter(predK, Level == "2"),
                  mapping = aes(Soil, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr),
                                group = Geno, col = Geno),
                  size = 1.2,
                  position = position_dodge(width = 0.5)) +
  xlab("Soil Treatment") +
  ylab("Potassium (mg/L)") +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

Kplot
ggsave(Kplot, filename = "figures/k.jpeg",
       height = 6, width =8, units = "in")

## Potassium Geno --------------------------------------------------------------
KpredG <- modpredG(KGeno)

KGplot <- ggplot() +
  geom_pointrange(data = KpredG,
                  mapping = aes(Geno, pred,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = KGeno$coefficients[1], col = "red") +
  xlab("Genotype") +
  ylab("log Potassium (mg/L)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))


ggsave(KGplot, filename = "figures/PotassiumGeno.jpeg",
       height = 6, width =8, units = "in")
## Magnesium Soil --------------------------------------------------------------
predMg <- modpredCat(Magnesium, Mg.boot)

# Create Plot

Mgplot <- ggplot() +
  geom_pointrange(data = filter(predMg, Level == "2"),
                  mapping = aes(Soil, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr),
                                group = Geno, col = Geno),
                  size = 1.2,
                  position = position_dodge(width = 0.5)) +
  xlab("Soil Treatment") +
  ylab("Magnesium (mg/L)") +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

Mgplot
ggsave(Mgplot, filename = "figures/Mg.jpeg",
       height = 6, width =8, units = "in")

## Magnesium Geno --------------------------------------------------------------
MgpredG <- modpredG(MgGeno)

MgGplot <- ggplot() +
  geom_pointrange(data = MgpredG,
                  mapping = aes(Geno, pred,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = MgGeno$coefficients[1], col = "red") +
  xlab("Genotype") +
  ylab("log Magnesium (mg/L)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))


ggsave(MgGplot, filename = "figures/MagnesiumGeno.jpeg",
       height = 6, width =8, units = "in")
## Sodium Soil -----------------------------------------------------------------
predNa <- modpredCat(Sodium, Na.boot)

# Create Plot

Naplot <- ggplot() +
  geom_pointrange(data = filter(predNa, Level == "2"),
                  mapping = aes(Soil, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr),
                                group = Geno, col = Geno),
                  size = 1.2,
                  position = position_dodge(width = 0.5)) +
  xlab("Soil Treatment") +
  ylab("Sodium (mg/L)") +
  scale_color_viridis_d() +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))
Naplot
ggsave(Naplot, filename = "figures/Na.jpeg",
       height = 6, width =8, units = "in")

## Sodium Geno -----------------------------------------------------------------
NapredG <- modpredG(NaGeno)

NaGplot <- ggplot() +
  geom_pointrange(data = NapredG,
                  mapping = aes(Geno, pred,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = NaGeno$coefficients[1], col = "red") +
  xlab("Genotype") +
  ylab("log Sodium (mg/L)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))


ggsave(NaGplot, filename = "figures/SodiumGeno.jpeg",
       height = 6, width =8, units = "in")
## Nitrate Soil ----------------------------------------------------------------
predNO3 <- modpred(Nitrate, NO3.boot)

# Create Plot

NO3plot <- ggplot() +
  geom_pointrange(data = filter(predNO3, Level == "1"),
                  mapping = aes(Soil, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  xlab("Soil Treatment") +
  ylab("Nitrate (mg/L)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

NO3plot
ggsave(NO3plot, filename = "figures/NO3.jpeg",
       height = 6, width =8, units = "in")

## Nitrate Geno ----------------------------------------------------------------
NO3predG <- modpredG(NO3Geno, NO3Geno.boot)

NO3Gplot <- ggplot() +
  geom_pointrange(data = NO3predG,
                  mapping = aes(Geno, preds,
                                ymin = lwr, ymax = upr),
                  size = 1.2,
                  position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = fixef(NO3Geno)[1], col = "red") +
  xlab("Genotype") +
  ylab("log Nitrate (mg/L)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

NO3Gplot

ggsave(NO3Gplot, filename = "figures/NitrateGeno.jpeg",
       height = 6, width =8, units = "in")
