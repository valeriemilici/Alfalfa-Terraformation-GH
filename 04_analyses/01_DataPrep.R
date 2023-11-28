### This script will prepare the alfalfa census data for performance-based
### analyses

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(lubridate) #manages dates
library(stringr) #deal with mismatched strings

#read in all performance data
dat <- read.csv("data/RawData/GCR_2A_CensusData.csv")
bmdat <- read.csv("data/RawData/biomassdata.csv")
licor <- read.csv("data/LiCorData/licordataclean.csv")
discharge <- read.csv("data/RawData/GCR_2A_DischargeData.csv")
anions <- read.csv("data/RawData/2A_anions.csv")
cations <- read.csv("data/RawData/2A_Discharge_Geochem_Master.csv")
TC <- read.csv("data/RawData/MP_Val_TC_Data_GCR_2023.csv")
### Create Performance Data ----------------------------------------------------

#prep licor data
licor$Photo <- ifelse(licor$Photo < 0 , 0, licor$Photo)
licor$Trmmol <- ifelse(licor$Trmmol < 0,0, licor$Trmmol)

licor2 <- licor %>% 
  group_by(pot) %>%
  summarise(meanPSN = mean(Photo),
            meanTPN = mean(Trmmol),
            leafarea = mean(leafarea)) %>%
  mutate(PSN = meanPSN/leafarea,
         TPN = meanTPN/leafarea)

#merge the files by pot code
dat2 <- left_join(dat, bmdat, by = "pot")
dat3 <- left_join(dat2, licor2, by = "pot")

#convert dates to date format
dat3$Date <- mdy(dat3$Date)
#Height must be numeric
dat3$Height <- as.numeric(dat3$Height)
dat3$Height.mm <- dat3$Height*100
#RGR function
RGRf <- function(h1, h2, t){
  h1_l = log(h1)
  h2_l = log(h2)
  RGR = (h2_l - h1_l)/t 
  return(RGR)
}

#create performance data
dat4 <- dat3 %>%
  #create new columns for model
  mutate(Soil = str_extract(pot, "\\w"),
         Geno = str_extract(pot, "\\w{2,}"),
         Level = substr(pot, 3,3),
         Location = str_c(table, rack, sep = "."),
         #probability of survival
         status = ifelse(Status == "D", 0,1)) %>% 
  #remove pots that don't have plants in them
  filter(Soil != "C" &
           Soil != "P" &
           Soil != "N") %>% 
  group_by(pot) %>%
  arrange(Date, .by_group = TRUE) %>%
  #add column for census interval
  mutate(censusint = as.numeric(Date - lag(Date, default = first(Date))),
         #add column for height previous census
         HtPrevCensus = lag(Height.mm),
         #compute RGR (cm/cm/day)
         RGR_d = RGRf(HtPrevCensus, Height.mm, censusint))

#create discharge experiment data ----------------------------------------------
dat5 <- left_join(discharge, dat3, by = "pot") %>%
  group_by(pot) %>%
  arrange(CensusNo) %>%
  filter(row_number() == n()) %>%
  mutate(Soil = str_extract(pot, "\\w"),
         Geno = str_extract(pot, "\\w{2,}"),
         Level = substr(pot, 3,3),
         Location = str_c(table, rack, sep = "."),
         #probability of survival
         status = ifelse(Status == "D", 0,1),
         et = potmass2 - potmass3,
         soil = case_when(Soil == "L" ~ "L",
                          Soil == "S" ~ "S",
                          Soil == "C" ~ "B",
                          Soil == "N" ~ "B"),
         plant = case_when(Soil == "L" ~ "P", 
                           Soil == "S" ~ "P",
                          TRUE ~ "B"),
         batch = case_when(dis_table == 2 ~ 1,
                           dis_table == 3 ~ 1,
                           dis_table == 4 ~ 1,
                           TRUE ~ 2)) %>%
  ungroup() %>%
  filter(et >= 0)

#Pore Water Geochemistry -------------------------------------------------------

#for some reason this is a character
anions$Cl.mg.L. <- as.numeric(anions$Cl.mg.L.)

anion <- anions %>% 
  mutate(Fluorine = 10* F.mg.L.,
         Chlorine = 10 * Cl.mg.L.,
         Nitrite = 10 * N02.mg.L.,
         Nitrate = 10*NO3.mg.L.,
         Phosphate = 10*PO4.mg.L.,
         Sulfate = 10* SO4.mg.L.) %>%
  group_by(pot) %>%
  summarise(Fluorine = mean(Fluorine),
            Chlorine = mean(Chlorine),
            Nitrite = mean(Nitrite),
            Nitrate = mean(Nitrate),
            Phosphate = mean(Phosphate),
            Sulfate = mean(Sulfate)) %>%
  mutate(Soil = str_extract(pot, "\\w"),
         Geno = str_extract(pot, "\\w{2,}"),
         Level = substr(pot, 3,3),
         soil = case_when(Soil == "L" ~ "L",
                          Soil == "S" ~ "S",
                          Soil == "C" ~ "B",
                          Soil == "N" ~ "B"),
         plant = case_when(Soil == "L" ~ "P", 
                           Soil == "S" ~ "P",
                           TRUE ~ "B"))
  #remove the RO water sample (it's got nothing in it)
  anion <- anion[-10,]
  
  #Cation, pH, and EC data
  cat <- cations %>% group_by(Sample.ID) %>%
    summarize(mu_Li = mean(Lithium),
              mu_Na = mean(Sodium),
              mu_K = mean(Potassium),
              mu_Mg = mean(Magnesium),
              mu_Ca = mean(Calcium),
              pH = pH,
              EC = EC..us.cm.) %>%
    rename(pot = Sample.ID) %>%
    mutate(Lithium = 10*mu_Li,
           Sodium = 10*mu_Na,
           Potassium = 10*mu_K,
           Magnesium = 10*mu_Mg,
           Calcium = 10*mu_Ca
           )
  
dischem <- merge(anion, cat, all = T, by = "pot")

# Fix TC data before merge
TC[23,1] <- "L-3-R24-03"
TC[52,1] <- "S-3-105-03"
TC <- TC[-138,]
TC[68,1] <- "L-3-G15-08"
TC[124,1] <- "N-3-NA-03"
TC[107,1] <- "L-3-INA-02"

TC$TC_ppm <- ifelse(TC$TC_ppm < 0, 0 , TC$TC_ppm)
# merge together
dischem2 <- merge(dischem, TC, all = T, by = "pot")
# add a column for cation batch because batch 2 is flawed
dischem2 <- dischem2 %>% mutate(cat_batch = case_when(Geno == "R24" ~ "2",
                                                Geno == "G74" ~ "2",
                                                Geno == "TAS" ~ "2",
                                                Geno == "YON" ~ "2",
                                                Geno == "INA" ~ "2",
                                                Geno == "K19" ~ "2",
                                                Geno == "CRI" ~ "2",
                                                Geno == "NA" ~ "2",
                                                TRUE ~ "1"))

dischem3 <- dischem2[,c(1, 8:12,26:27,2:7,18:25)]

#add Dmass_tot to the spreadsheet
dmass <- dat5[,c(3,11)]

dischem4 <- left_join(dischem3, dmass, by = "pot")

#remove L-3-K19-03 (sample does not exist)
dischem4 <- dischem4[-58,]

#Generate Metadata File --------------------------------------------------------

metadata<- dat5 %>% select(pot, Soil, Level, Geno, Location)

#Export the finished dataframes ------------------------------------------------
write.csv(dat4, file = "data/performance.csv")
write.csv(dat5, file = "data/discharge.csv")
write.csv(dischem4, file = "data/2A_porewater_geochemistry.csv")
write.csv(metadata, file = "data/metadata.csv")
