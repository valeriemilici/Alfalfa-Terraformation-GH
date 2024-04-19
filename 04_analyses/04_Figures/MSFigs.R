### The figures used in the Manuscript

### Initialize workspace--------------------------------------------------------
rm(list =ls())
library(lme4)
library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(patchwork) #joins plots nicely
library(ggeffects) #plot directly from model without bootstrap
library(stats)
library(parallel)

### Figure 1 (Microbe Effects on Plant Performance) ----------------------------

#the models to plot from:
mod1 <- readRDS("04_analyses/02_models/Output/RGR.RDS")
mod2 <- readRDS("04_analyses/02_models/Output/Biomass.RDS")
mod3 <- readRDS("04_analyses/02_models/Output/Surv.RDS")

# Panel A

RGRpreds <- ggemmeans(mod1, c("Level", "Soil"), type = "fixed")

Fig1A <- ggplot() +
  geom_pointrange(data = RGRpreds,
                  mapping = aes(x, predicted,
                                ymin = conf.low,
                                ymax = conf.high,
                                group = group,
                                col = group),
                  position = position_dodge(width = 0.2)) +
  xlab("") +
  ylab("Relative Growth Rate \n(mm/mm/day)") +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Fig1A

# Panel B

BMpreds <- ggemmeans(mod2, c("Level", "Soil"), type = "fixed")

Fig1B <- ggplot() +
  geom_pointrange(data = BMpreds,
                  mapping = aes(x, predicted,
                                ymin = conf.low,
                                ymax = conf.high,
                                group = group,
                                col = group),
                  position = position_dodge(width = 0.2)) +
  xlab("") +
  ylab("Total Dry Biomass (g)") +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Fig1B

# Panel C

Survpreds <- ggemmeans(mod3, c("Level", "Soil"), type = "fixed")

Fig1C <- ggplot() +
  geom_pointrange(data = Survpreds,
                  mapping = aes(x, predicted,
                                ymin = conf.low,
                                ymax = conf.high,
                                group = group,
                                col = group),
                  position = position_dodge(width = 0.2)) +
  xlab("Percent of Water Holding Capacity") +
  ylab("Probability of Survival") +
  scale_x_discrete(labels = c("1" = "30-50", "2" = "55-75", "3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c("L" = "Live Soil", "S" = "Sterile Soil"),
                        begin = 0.2,
                        end = 0.8) +
  theme_classic(12)

Fig1C

# Combine Panels

Fig1 <- Fig1A/Fig1B/Fig1C +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

Fig1

# Save Figure
ggsave(plot = Fig1, filename = "04_analyses/04_Figures/ManuscriptFigs/Figure1.png",
       width = 4, height = 6, units = "in")

### Figure 2 (Root Mass Fraction/ Biomass Proportion) --------------------------

mod4 <- readRDS("04_analyses/02_models/Output/RMF.RDS")

RMFpreds <- ggemmeans(mod4, "Soil", type = "fixed")

RMFpreds$shoot <- 1 - RMFpreds$predicted

RMFpreds2 <- data.frame(RMFpreds) %>% dplyr::select(x, predicted, shoot)

RMFpreds2 <- pivot_longer(RMFpreds2, cols = c("predicted", "shoot"))
RMFpreds2$name <- factor(RMFpreds2$name, levels = c("shoot", "predicted"))


Fig2 <- ggplot() +
  geom_bar(RMFpreds2, mapping = aes(fill = name, y = value, x = x),
           position = "stack", stat = "identity") +
  geom_errorbar(RMFpreds, mapping = aes(x = x,
                                        ymin = conf.low, ymax = conf.high),
                width = 0.1) +
  xlab("Treatment") +
  ylab("Biomass Proportion") +
  scale_x_discrete(labels = c("L" = "Live Soil", "S" = "Sterile Soil")) +
  scale_fill_viridis_d("",
                       labels = c("shoot" = "Shoot Mass", "predicted" = "Root Mass"),
                        begin = 0.2,
                        end = 0.8,
                       direction = -1) +
  theme_classic(12)

Fig2

ggsave(plot = Fig2, filename = "04_analyses/04_Figures/ManuscriptFigs/Figure2.png")


### Figure 3 (Genotype Specific Relationship with Soil Microbes and Biomass)----

mod5 <- readRDS("04_analyses/02_models/Output/BiomassGS.RDS")

genoin1 <- c("CRI-CRI", "G15-G15", "K19-K19", "MOC-MOC", "TAS-TAS", "VIR-VIR",
            "YON-YON")
m5preds <- data.frame(test_predictions(mod5, c("Geno", "Soil"))) %>%
  filter(Geno %in% genoin)

Fig3 <- ggplot(m5preds, aes(Geno, (-1* Contrast/10))) +
  geom_pointrange(mapping = aes(ymin = (-1* conf.low/10), ymax = (-1* conf.high/10))) + 
  labs(x = "Population",
       y = "Microbe Effect on Biomass") +
  scale_x_discrete(labels = c("CRI-CRI" = "1",
                              "G15-G15" = "2",
                              "K19-K19" = "3",
                              "MOC-MOC" = "4",
                              "TAS-TAS" = "5",
                              "VIR-VIR" = "6",
                              "YON-YON" = "7")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(12) 

Fig3

ggsave(plot = Fig3, filename = "04_analyses/04_Figures/ManuscriptFigs/Figure3.png",
       height = 4, width = 6, units = "in")

### Figure 4 (EC and pH) -------------------------------------------------------

mod6 <- readRDS("04_analyses/02_models/Output/ECSoil.RDS")
mod7 <- readRDS("04_analyses/02_models/Output/pHSoil.RDS")

#Panel A

pHpreds <- ggemmeans(mod7, terms = "Soil", type = "fixed")

pHpreds$x <- factor(pHpreds$x, c("C", "N", "S", "L"))

Fig4A <- ggplot() +
  geom_pointrange(data = pHpreds,
                  mapping = aes(x, predicted,
                                ymin = conf.low, ymax = conf.high),
                  size = 1,
                  position = position_dodge(width = 0.2)) +
  xlab("") +
  ylab("pH") +
  theme_classic(12) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

Fig4A

#Panel B

ECpreds <- ggemmeans(mod6, terms = "Soil", type = "fixed")

ECpreds$x <- factor(ECpreds$x, c("C", "N", "S", "L"))

Fig4B <- ggplot() +
  geom_pointrange(data = ECpreds,
                  mapping = aes(x, predicted,
                                ymin = conf.low, ymax = conf.high),
                  size = 1,
                  position = position_dodge(width = 0.2)) +
  xlab("Soil Biota Treatment") +
  ylab("Electrical Conductivity (EC) ") +
  scale_x_discrete(labels = c("C" = "Abiotic \n Control", "N" = "Microbes \nOnly",
                              "S" = "Plant \nOnly", "L" = "Plant & \n Microbes")) +
  theme_classic(12)

Fig4B

# Figure 4

Fig4 <- Fig4A/Fig4B + plot_annotation(tag_levels = "A")
Fig4

ggsave(plot = Fig4, filename = "04_analyses/04_Figures/ManuscriptFigs/Figure4.png")

### Figure 5 (Broad PW Geochem Trends) -----------------------------------------
# load all of those geochemistry models
mod8 <- readRDS("04_analyses/02_models/Output/ICSoil.RDS")
mod9 <- readRDS("04_analyses/02_models/Output/OCSoil.RDS")
mod10 <- readRDS("04_analyses/02_models/Output/NO3Soil.RDS")
mod11 <- readRDS("04_analyses/02_models/Output/NH4Soil.RDS")
mod12 <- readRDS("04_analyses/02_models/Output/NaSoil.RDS")
mod13 <- readRDS("04_analyses/02_models/Output/CaSoil.RDS")
mod14 <- readRDS("04_analyses/02_models/Output/LiSoil.RDS")
mod15 <- readRDS("04_analyses/02_models/Output/MgSoil.RDS")
mod16 <- readRDS("04_analyses/02_models/Output/KSoil.RDS")
mod17 <- readRDS("04_analyses/02_models/Output/FSoil.RDS")
mod18 <- readRDS("04_analyses/02_models/Output/BrSoil.RDS")
mod19 <- readRDS("04_analyses/02_models/Output/ClSoil.RDS")
mod20 <- readRDS("04_analyses/02_models/Output/SO4Soil.RDS")
#mod21 <- readRDS("04_analyses/02_models/Output/PO4Soil.RDS")

#Note- must remove PO4 model because some categories lack data and
# can't be used for predictions

# A function that generate model predictions for each treatment,
# computes the difference relative to the abiotic control treatment,
# and affixes the p-value for each comparison.
heatpred <- function(mod){
  
  modpred <- data.frame(ggemmeans(mod, c("Soil", "Level"))) %>%
  dplyr::select(x, predicted, group) %>%
  rename(Soil = x, Level = group)

control <- modpred %>% filter(Soil == "C")

modpred <- modpred %>% filter(Soil != "C")

modpred <- left_join(modpred, control, by = "Level") %>%
  group_by(Level) %>%
  mutate(prop.diff = ((predicted.x - predicted.y) / predicted.y)) %>%
  dplyr::select(Soil.x, Level, prop.diff) %>%
  rename(Soil = Soil.x)

levelin <- c("1-1", "2-2", "3-3")
soilin <- c("C-S", "C-L", "C-N")

pv <- data.frame(test_predictions(mod, terms = c("Soil", "Level"))) %>%
filter(Level %in% levelin & Soil %in% soilin) %>%
  mutate(Soil2 = str_extract(Soil, '\\w$'),
         Level2 = str_extract(Level, '\\w$')) %>%
  dplyr::select(Soil2, Level2, p.value) %>%
  rename(Soil = Soil2, Level = Level2)

pred <- left_join(modpred, pv, by = c("Soil", "Level")) 

return(pred)
}

# create a list of the models to run through the function
modlist <- list(mod8, mod9, mod10, mod11, mod12,
                mod13, mod14, mod15, mod16, mod17,
                mod18, mod19, mod20
             #   , mod21
                )

# run list of models through the function and
# transform output so that it is usable downstream
predout <- data.frame(apply(sapply(modlist, heatpred), 1, unlist)) %>%
  mutate(Solute = rep(c("Inorganic Carbon",
                        "Organic Carbon",
                        "Nitrate",
                        "Ammonium",
                        "Sodium",
                        "Calcium",
                        "Lithium",
                        "Magnesium",
                        "Potassium",
                        "Fluoride",
                        "Bromide",
                        "Chloride",
                        "Sulfate"
                       # ,"Phosphate"
                        ), each = 9))

#for some reason this bit messed up Lithium when it was run as part of
# the previous chunk. It's totally fine when separated. Weird. 
predout2 <- predout %>% mutate( prop.diff.t = ifelse(p.value > 0.1, 0, prop.diff),
         p.value.t = case_when(p.value < 0.001 ~ "***",
                               p.value < 0.01 ~ "**",
                               p.value < 0.05 ~ "*",
                               p.value < 0.1 ~ ".",
                               .default = " "))
# these columns must be numeric
predout2$prop.diff <- as.numeric(predout2$prop.diff)
predout2$p.value <- as.numeric(predout2$p.value)

# Make the figure

predout2$Soil <- factor(predout2$Soil, c("N", "S", "L"))

predout2$Solute <- factor(predout2$Solute, c("Inorganic Carbon",
                                             "Organic Carbon",
                                             "Nitrate",
                                             "Ammonium",
                                             "Sodium",
                                             "Calcium",
                                             "Lithium",
                                             "Magnesium",
                                             "Potassium",
                                             "Fluoride",
                                             "Bromide",
                                             "Chloride",
                                             "Sulfate"))
Fig5 <- ggplot(predout2, aes(x = Soil, y = Solute, fill = log1p(prop.diff))) +
  geom_tile(color = "black") +
  facet_wrap(~Level, labeller = as_labeller(c("1" = "30-50% of WHC",
                                              "2" = "55-75% of WHC",
                                              "3" = "80-100% of WHC"))) +
    geom_text(aes(label = p.value.t),
              size = 8) +
  xlab("Soil Biota Treatment") +
    scale_fill_gradient2("",
                         low = "blue", 
                        mid = "white",
                        high = "red") +
  scale_x_discrete(labels = c("N" = "Microbes \nOnly",
                              "S" = "Plant \nOnly",
                              "L" = "Plant & \n Microbes")) +
  coord_fixed() +
  theme_classic(12)

Fig5

ggsave(plot = Fig5,
       filename = "04_analyses/04_Figures/ManuscriptFigs/Figure5.png")


### Figure 6 (Specific Highlighted Trends) -------------------------------------

# Panel A- Sodium

levelin <- c("1-1", "2-2", "3-3")
soilin <- c("C-N", "C-L", "C-S")

pvNa <- data.frame(test_predictions(mod12, terms = c("Soil", "Level"))) %>%
  filter(Soil %in% soilin &
           Level %in% levelin)

pvNa$Soil <- factor(pvNa$Soil, levels = c("C-N", "C-S", "C-L"))

panel6A <- ggplot(data = pvNa,
                  mapping = aes(Level, (Contrast* -1),
                                group = Soil, color = Soil)) +
  geom_line() +
  geom_pointrange(mapping = aes(ymin = (conf.high * -1),
                                ymax = (conf.low * -1)),
                  position = position_dodge(width = 0.1)) +
  geom_hline(yintercept = 0,
             linetype = "dotted") +
  
  xlab("") +
  ylab("Sodium (mg/L)") +
  scale_x_discrete(labels = c("1-1" = "30-50",
                              "2-2" = "55-75",
                              "3-3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c(
                          "C-N" = "Microbes Only",
                          "C-S" = "Plant Only",
                          "C-L" = "Plant & \nMicrobes"),
                        begin = 0.2) +
  theme_classic(12)

panel6A

#Panel 6B - Lithium

pvLi <- data.frame(test_predictions(mod14, terms = c("Soil", "Level"))) %>%
  filter(Soil %in% soilin &
           Level %in% levelin)

pvLi$Soil <- factor(pvLi$Soil, levels = c("C-N", "C-S", "C-L"))

panel6B <- ggplot(data = pvLi,
                  mapping = aes(Level, (Contrast* -1),
                                group = Soil, color = Soil)) +
  geom_line() +
  geom_pointrange(mapping = aes(ymin = (conf.high * -1),
                                ymax = (conf.low * -1)),
                  position = position_dodge(width = 0.1)) +
  geom_hline(yintercept = 0,
             linetype = "dotted") +
  
  xlab("") +
  ylab("Lithium (mg/L)") +
  scale_x_discrete(labels = c("1-1" = "30-50",
                              "2-2" = "55-75",
                              "3-3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c(
                          "C-N" = "Microbes Only",
                          "C-S" = "Plant Only",
                          "C-L" = "Plant & \nMicrobes"),
                        begin = 0.2) +
  theme_classic(12)


panel6B
#Figure Panel 6C - Magnesium
                                      
pvMg <- data.frame(test_predictions(mod15, terms = c("Soil", "Level"))) %>%
  filter(Soil %in% soilin &
           Level %in% levelin)

pvMg$Soil <- factor(pvMg$Soil, levels = c("C-N", "C-S", "C-L"))

panel6C <- ggplot(data = pvMg,
                  mapping = aes(Level, (Contrast * -1),
                                group = Soil, color = Soil)) +
  geom_line() +
  geom_pointrange(mapping = aes(ymin = (conf.high * -1),
                                ymax = (conf.low * -1)),
                  position = position_dodge(width = 0.1)) +
  geom_hline(yintercept = 0,
             linetype = "dotted") +
 
  xlab("") +
  ylab("Magnesium (mg/L)") +
  scale_x_discrete(labels = c("1-1" = "30-50",
                              "2-2" = "55-75",
                              "3-3" = "80-100")) +
  scale_color_viridis_d("",
                        labels = c(
                                   "C-N" = "Microbes Only",
                                   "C-S" = "Plant Only",
                                   "C-L" = "Plant & \nMicrobes"),
                        begin = 0.2) +
  theme_classic(12)

panel6C

# Join Together

Figure6 <- panel6A + panel6B + panel6C +
  plot_annotation(tag_levels = "A",
                  caption = "Percent of Total Water Holding Capacity",
                  theme = theme(plot.caption = element_text(hjust = 0.45, size = 12))) +
  plot_layout(guides = "collect")

ggsave(plot = Figure6, filename = "04_analyses/04_Figures/ManuscriptFigs/Figure6.png",
       width = 10, height = 3, units = "in")

### Figure NA (Solutes responsive to moisture only) ------------------------------

#Note- there is no solute that is responsive to soil moisture alone

### Figure 7 (GS microbe effects on IC) ----------------------------------------
mod22 <- readRDS("04_analyses/02_models/Output/ICGeno.RDS")

genoin <- c("CRI-CRI", "G15-G15", "K19-K19", "MOC-MOC", "TAS-TAS", "VIR-VIR",
            "YON-YON")
m22preds <- data.frame(test_predictions(mod22, c("Geno", "Soil"))) %>%
  filter(Geno %in% genoin)

Fig7 <- ggplot(m22preds, aes(Geno, Contrast)) +
  geom_pointrange(mapping = aes(ymin = conf.low, ymax = conf.high)) + 
  labs(x = "Population",
       y = "Microbe Effect on \nInorganic Carbon (ppm)") +
  scale_x_discrete(labels = c("CRI-CRI" = "1",
                              "G15-G15" = "2",
                              "K19-K19" = "3",
                              "MOC-MOC" = "4",
                              "TAS-TAS" = "5",
                              "VIR-VIR" = "6",
                              "YON-YON" = "7")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic(12) 

Fig7
ggsave(plot = Fig7, filename = "04_analyses/04_Figures/ManuscriptFigs/Figure7.png",
       height = 4, width = 6, units = "in")

### Figure NA (GS water effects) ----------------------------------------

#all of these are terrible and don't add to storytelling

### Supplemental Table 1 ------------------------

tablefun<- function(mod){
levelin <- c("1-1", "2-2", "3-3")
soilin <- c("C-S", "C-L", "C-N")
pv <- data.frame(test_predictions(mod, terms = c("Soil", "Level"))) %>%
  filter(Level %in% levelin & Soil %in% soilin) %>%
  mutate(Soil2 = str_extract(Soil, '\\w$'),
         Level2 = str_extract(Level, '\\w$')) %>%
  dplyr::select(!c(Soil, Level)) %>%
  rename(Soil = Soil2, Level = Level2, Est.Diff = Contrast)

return(pv)

}
modlist <- list(mod8, mod9, mod10, mod11, mod12,
                mod13, mod14, mod15, mod16, mod17,
                mod18, mod19, mod20
                #   , mod21
)

# run list of models through the function and
# transform output so that it is usable downstream
tableout <- data.frame(apply(sapply(modlist, tablefun), 1, unlist)) %>%
  mutate(Solute = rep(c("Inorganic Carbon",
                        "Organic Carbon",
                        "Nitrate",
                        "Ammonium",
                        "Sodium",
                        "Calcium",
                        "Lithium",
                        "Magnesium",
                        "Potassium",
                        "Fluoride",
                        "Bromide",
                        "Chloride",
                        "Sulfate"
                        # ,"Phosphate"
  ), each = 9))

write.csv(tableout, file = "04_analyses/04_Figures/ManuscriptFigs/supptable1.csv")
