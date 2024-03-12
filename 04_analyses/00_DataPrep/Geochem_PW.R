### Collates all porewater geochemistry files
### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(lubridate) #manages dates
library(stringr) #deal with mismatched strings

#Types of files to collate include: 1. Cations, 2. Anions, 3. TC/OC/IC, 4. pH/EC
# Cation
b1 <- read.csv("data/00_PrepData/MP_Val_Cation_Reruns_Sept_7_2023.csv")
b2<- read.csv("data/00_PrepData/MP_Val_Cation_Reruns_Sept_11_2023.csv")
b3<- read.csv("data/00_PrepData/MP_Val_Cation_Reruns_Sept_18_2023.csv")
b4<- read.csv("data/00_PrepData/MP_Val_Cation_Reruns_Sept_25_2023.csv")
b1Cats <- read.csv("data/00_PrepData/Partial_Cation_Clean.csv")
# Anion
anion <- read.csv("data/00_PrepData/2A_anions.csv")
anion_batch <- read.csv("data/00_PrepData/AnionBatch.csv")
#Carbon
TC <- read.csv("data/00_PrepData/TC_dat.csv")
NPOC <- read.csv("data/00_PrepData/NPOC_Dat.csv")
#pH and EC
EC <- read.csv("data/00_PrepData/2A_discharge_PhEC.csv")
#Step 1: Cations ---------------------------------------------------------------
#merge together
cation <- rbind(b1,b2,b3,b4)

#create "rep" column
cation$rep <- str_extract(cation$pot, '\\d$')

#remove "rep" value from pot code
cation$pot <- str_sub(cation$pot, end = -3)

#condense reps into one mean value
CatDat <- cation %>% group_by(pot, batch) %>%
  summarise(mu_Li = mean(lithium, na.rm = T),
            mu_Na = mean(sodium, na.rm = T),
            mu_NH4 = mean(ammonium, na.rm = T),
            mu_K = mean(potassium, na.rm = T),
            mu_Mg = mean(magnesium, na.rm = T),
            mu_Ca = mean(calcium, na.rm = T))
#Make identical to CatDat
b1Cat <- b1Cats %>%
  rename(pot = Sample) %>%
  mutate(batch = as.integer(5),
         mu_NH4 = as.numeric(NA)) %>%
  dplyr::select(c(1,7,2,3,8,4:6))

#object to be used with the final merge
Cation <- rbind(CatDat, b1Cat)

Cation <- rename(Cation, batch_cat = batch)
Cation[137,1] <- "S-3-G15-02"

# Step 2: Anions ---------------------------------------------------------------
# Compute mean anion
# use this for merging
# multiplying by 10 here to account for dilution. Confirm with geochemists that
# the scale of values seems appropriate.
Anion <- anion %>%
  group_by(pot) %>%
  summarize(mu_F = 10*mean(F.mg.L., na.rm = T),
            mu_Cl = 10*mean(as.numeric(Cl.mg.L.), na.rm = T),
            mu_NO2 = 10*mean(N02.mg.L., na.rm = T),
            mu_Br = 10*mean(Br.mg.L., na.rm = T),
            mu_NO3 = 10*mean(NO3.mg.L., na.rm = T),
            mu_PO4 = 10*mean(PO4.mg.L., na.rm = T),
            mu_SO4 = 10*mean(SO4.mg.L., na.rm = T)) 
# combine batch information with the data

Anion2 <- merge(Anion, anion_batch, by = "pot", all = T)

# manually add in anion information to row 14 (L-1-g74-01)
# data/Extras/ExcelFiles/Anions/Anions/*April_13_2023

Anion2[14,2] <- 5.767 #fluoride
Anion2[14,3] <- 139.288 #chloride
Anion2[14,5] <- 3.961 #bromide
Anion2[14,6] <- 10.217 #nitrate
Anion2[14,8] <- 208.046 #sulfate


# Step 3: TC/OC/IC -------------------------------------------------------------

Carbon <- merge(TC,NPOC, by = "pot", all = T)

Carbon <- Carbon %>%
  mutate(IC = TC_ppm - NPOC,
         batch_IC = paste(batch_TC, batch_NPOC, sep = ".")) 

Carbon$IC <- ifelse(Carbon$IC < 0, 0, Carbon$IC)

Carbon[61,1] <- "L-3-INA-02"
#Step 4: EC/pH -----------------------------------------------------------------

EC <- EC %>%
  dplyr::select(2,4,5) %>%
  rename(EC = EC..us.cm.)

#Step 5: Merge them all together -----------------------------------------------
PW_1 <- merge(Cation, Anion2, all = T)
PW_2 <- merge(PW_1, Carbon, all = T)
PW_3 <- merge(PW_2, EC, all = T)

#nonsense rows that are empty and the samples don't exist
outs <- c("L-2-CRI-02", "L-2-CRI-07", "L-3-R25-03", "GH RO SAMPLE")

PW_3 <- filter(PW_3, !pot %in% outs)
#save
write.csv(PW_3,"data/01_PrepData/PW_Geochemistry.csv")
