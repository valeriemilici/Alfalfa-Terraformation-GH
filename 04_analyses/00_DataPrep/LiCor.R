### This code will take the individual LiCor data sheets and merge them into
### a master data sheet that can be used for analysis.

### Intialize Worksheet --------------------------------------------------------
rm(list =ls())
#packages
library(tidyverse) #always useful

#data
leafarea <- read.csv("data/00_PrepData/LiCorData/LiCor_LeafArea.csv")
dat1<- read.csv("data/00_PrepData/LiCorData/09_26_22_licor.csv")
dat2<- read.csv("data/00_PrepData/LiCorData/09_28_22_licor.csv")
dat3<- read.csv("data/00_PrepData/LiCorData/09_30_22_licor.csv")
dat4<- read.csv("data/00_PrepData/LiCorData/10_03_22_licor.csv")
dat5<- read.csv("data/00_PrepData/LiCorData/10_05_22_licor.csv")
dat6<- read.csv("data/00_PrepData/LiCorData/10_05_22_licor_1.csv")
dat7<- read.csv("data/00_PrepData/LiCorData/10_07_22_licor.csv")
dat8<- read.csv("data/00_PrepData/LiCorData/10_09_22_licor.csv")

### Merge LiCor Obs ------------------------------------------------------------

LiCorDat<- rbind(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8)
#View(LiCorDat)
#looks good!

### Add in leaf area -----------------------------------------------------------
leafarea2 <- leafarea %>% dplyr::select(1:2) %>%
  rename(leafarea = area.mm.sq.)
leafarea2$pot <- toupper(leafarea2$pot) #just in case

leafarea2[124,1] <- "S-1-YON-03" #improperly named

# Convert negative values to zero
LiCorDat$Photo <- ifelse(LiCorDat$Photo < 0, 0, LiCorDat$Photo)
LiCorDat$Trmmol <- ifelse(LiCorDat$Trmmol < 0, 0, LiCorDat$Trmmol)
# Looks like some pot names were lowercase (will create downstream issues)
LiCorDat$pot <- toupper(LiCorDat$pot)

#Simplify LiCorDat to better diagnose merge issues
LiCorDat2 <- LiCorDat %>%
  group_by(pot) %>%
  summarize(mean_PSN = mean(Photo),
            mean_TPN = mean(Trmmol))

LiCorDat2[13,1] <- "L-2-515-05"
#combine the datasets
FullDat <- merge(LiCorDat2, leafarea2, by = "pot")

FullDat2 <- FullDat %>%
  mutate(std_PSN = mean_PSN/leafarea,
         std_TPN = mean_TPN/leafarea)

### save file ------------------------------------------------------------------
write.csv(FullDat2, file = "data/01_PrepData/LiCor.csv")
