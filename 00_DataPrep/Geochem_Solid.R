### Collates all solid geochemistry files
### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(lubridate) #manages dates
library(stringr) #deal with mismatched strings

TC <- read.csv("data/00_PrepData/SolidSamples_C_N_PriorityOnly.csv")
TOC <- read.csv("data/00_PrepData/solid_TOC.csv")

# combine datasets and then calculate IC

TC1 <- TC %>%
  rename(pot = Sample.ID) %>%
  group_by(pot, Acid.Treatment) %>%
  summarize(mu_TC = mean(TC_ug.mg, na.rm = T),
            mu_TN = mean(TN_ug.mg, na.rm = T),
            batch_STC = mean(batch))

TC1$Acid.Treatment <- toupper(TC1$Acid.Treatment)

#TC is non-acid treatment only
TC2 <- TC1 %>% 
  filter(Acid.Treatment == "NO") %>%
  ungroup() %>%
  select(1, 3:5)

#clean up and make TOC data sets identical
#some OC data were in the other dataset, filter and rename
TOC2 <- TC1 %>%
  filter(Acid.Treatment == "YES") %>%
  ungroup() %>%
  select(pot, mu_TC, batch_STC) %>%
  rename(mu_TOC = mu_TC,
         batch_STOC = batch_STC)

TOC1 <- TOC %>%
  rename(pot = ID) %>%
  group_by(pot) %>%
  summarize(mu_TOC = mean(ug_OC_mg),
            batch_STOC = mean(batch_SOC)) %>%
  filter(batch_STOC > 0)

TOC3 <- rbind(TOC1, TOC2) #create unified TOC dataset/column

#combine data sets
Solid_Carbon <- merge(TC2, TOC3, by = "pot") #looks good!

#calculate IC

Carbon_s <- Solid_Carbon %>%
  mutate(IC_s = mu_TC - mu_TOC,
         batch_SIC = paste(batch_STC, batch_STOC, sep = ".")) %>%
  rename(TC_s = mu_TC,
         TOC_s = mu_TOC,
         TN_s = mu_TN) %>%
  select(1:3,5,7,4,6,8)

#negative carbon isn't possible. Convert to zero
Carbon_s$IC_s <- ifelse(Carbon_s$IC_s <= 0, 0, Carbon_s$IC_s) 

#save file
write.csv(Carbon_s, file = "data/01_PrepData/Solid_Geochemistry_Priority.csv")
