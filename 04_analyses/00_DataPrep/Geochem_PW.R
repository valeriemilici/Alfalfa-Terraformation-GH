### Collates all porewater geochemistry files
### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(lubridate) #manages dates
library(stringr) #deal with mismatched strings

#Types of files to collate include: 1. Cations, 2. Anions, 3. TC/OC/IC, 4. pH/EC
# Cation

cation <- read.csv("data/00_PrepData/Cation_Clean_Complete.csv")
cat_charge <- read.csv("data/00_PrepData/totalcatconcentration.csv")

# Anion

anion <- read.csv("data/00_PrepData/Anions_Clean_Complete.csv")

#delete this data or move to new folder
#anion <- read.csv("data/00_PrepData/2A_anions.csv")
#anion_batch <- read.csv("data/00_PrepData/AnionBatch.csv")

#Carbon
TC <- read.csv("data/00_PrepData/TC_dat.csv")
NPOC <- read.csv("data/00_PrepData/NPOC_Dat.csv")
#Nitrogen
TN <- read.csv("data/00_PrepData/TN_PW.csv")
#pH and EC
EC <- read.csv("data/00_PrepData/2A_discharge_PhEC.csv")
#Step 1: Cations ---------------------------------------------------------------

# add concentration to other cat data
cation <- left_join(cation, cat_charge, by = "pot")

#create "rep" column
cation$rep <- str_extract(cation$pot, '\\d$')

#remove "rep" value from pot code
cation$pot <- str_sub(cation$pot, end = -3)

#split data into two sets of batches

cat2<- cation %>% filter(batch_cation >4)
cation <- cation %>% filter(batch_cation <= 4)


#condense reps into one mean value
CatDat1 <- cation %>% group_by(pot, batch_cation) %>%
  summarise(mu_Li = mean(Li_diluted, na.rm = T)*10,
            mu_Na = mean(Na_diluted, na.rm = T)*10,
            mu_NH4 = mean(NH4_diluted, na.rm = T)*10,
            mu_K = mean(K_diluted, na.rm = T)*10,
            mu_Mg = mean(Mg_diluted, na.rm = T)*10,
            mu_Ca = mean(Ca_diluted, na.rm = T)*10,
            mu_charge = mean(Cat_Charge, na.rm = T)) #charge pre-adjusted

CatDat2 <- cat2 %>% group_by(pot, batch_cation) %>%
 summarise(mu_Li = mean(Li_diluted, na.rm = T)*10*2.6,
           mu_Na = mean(Na_diluted, na.rm = T)*10*2.6,
           mu_NH4 = mean(NH4_diluted, na.rm = T)*10*2.6,
           mu_K = mean(K_diluted, na.rm = T)*10*2.6,
           mu_Mg = mean(Mg_diluted, na.rm = T)*10*2.6,
            mu_Ca = mean(Ca_diluted, na.rm = T)*10*2.6,
           mu_charge= mean(Cat_Charge, na.rm =T)*2.6)

#object to be used with the final merge
 Cation <- rbind(CatDat1, CatDat2)

# Step 2: Anions ---------------------------------------------------------------
# Compute mean anion
# use this for merging
# multiplying by 10 here to account for dilution.

#create rep column and remove rep from pot code, as above
anion$rep <- str_extract(anion$pot, '\\d$')

anion$pot <- str_sub(anion$pot, end = -3)

Anion <- anion %>%
  group_by(pot, batch_anion) %>%
  summarize(mu_F = 10*mean(F_diluted, na.rm = T),
            mu_Cl = 10*mean(Cl_diluted, na.rm = T),
            mu_NO2 = 10*mean(nitrite_diluted, na.rm = T),
            mu_Br = 10*mean(br_diluted, na.rm = T),
            mu_NO3 = 10*mean(no3_diluted, na.rm = T),
            mu_PO4 = 10*mean(po4_diluted, na.rm = T),
            mu_SO4 = 10*mean(so4_diluted, na.rm = T)) 


# fix sample name
Anion[14,1] <- "L-1-G74-01"

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

#Step 5: TN --------------------------------------------------------------------
TN <- TN %>% rename(TN.ppm = TN..ppm.)
#Step 5: Merge them all together -----------------------------------------------
PW_1 <- merge(Cation, Anion, all = T)
PW_2 <- merge(PW_1, Carbon, all = T)
PW_3 <- merge(PW_2, EC, all = T)
PW_4 <- merge(PW_3, TN, all = T)

#nonsense rows that are empty and the samples don't exist
outs <- c("L-2-CRI-02", "L-2-CRI-07", "L-3-R25-03", "GH RO SAMPLE", "GH-RO")

PW_4 <- filter(PW_4, !pot %in% outs)
#save
write.csv(PW_4,"data/01_PrepData/PW_Geochemistry.csv")
