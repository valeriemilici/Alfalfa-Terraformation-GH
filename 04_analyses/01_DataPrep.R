### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(lubridate) #manages dates
library(stringr) #deal with mismatched strings

#read in all performance data

#a time series data set with height and survival over duration of exp.
#this is the "lead" .csv and other files will be merged within it.
census <- read.csv("data/01_PrepData/GCR_2A_CensusData.csv")

#the above and belowground biomass for all plants that survived
bm <- read.csv("data/01_PrepData/biomassdata.csv")

#data from the discharge experiment, can be used to calculate ET
discharge <- read.csv("data/01_PrepData/GCR_2A_DischargeData.csv")

#Photosynthetic and transpiration rates
licor <- read.csv("data/01_PrepData/LiCor.csv")

#Collated porewater geochemistry data (complete)
pw_geo <- read.csv("data/01_PrepData/PW_Geochemistry.csv")

#Carbon and Nitrogen information in the soil (incomplete)
s_geo <- read.csv("data/01_PrepData/Solid_Geochemistry_Priority.csv")

#root length of the alfalfa
root <- read.csv("data/01_PrepData/RootLength.csv")
### Create Performance Data ----------------------------------------------------

#We will create two different spreadsheets: 
# 1. clean up "census" to for a model-ready time series.
# 2. merge all other data sets together to create a spreadsheet of all of the
# final measurements on the surviving plants. 

# 1. Creating a clean time series from "census" --------------------------------

#convert dates to date format
census$Date <- mdy(census$Date)
#Height must be numeric
census$Height <- as.numeric(census$Height)
# mm is better scale for model estimates
census$Height.mm <- census$Height*100
#RGR function
RGRf <- function(h1, h2, t){
  h1_l = log(h1)
  h2_l = log(h2)
  RGR = (h2_l - h1_l)/t 
  return(RGR)
}

#create performance data
census1 <- census %>%
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
         #compute RGR (mm/mm/day)
         RGR_d = RGRf(HtPrevCensus, Height.mm, censusint)) %>%
  ungroup() %>%
  dplyr::select(1,2,3,13:15,17,6,11,5,12,18:20,7,8,16) 
#Looks good!

# 2. Merge together all final measurements -------------------------------------

#trim down and clean up full census data
dat <- census %>%
  mutate(Soil = str_extract(pot, "\\w"),
         Geno = str_extract(pot, "\\w{2,}"),
         Level = substr(pot, 3,3),
         Location = str_c(table, rack, sep = "."),
         #probability of survival
         status = ifelse(Status == "D", 0,1)) %>%
  #only need final census data
  filter(CensusNo == 6) %>%
  dplyr::select(3,17,5,12,11,6:8,13:16)

#discharge data can only tell us ET (following discussion with Hannes B.)
dis1 <- discharge %>%
  mutate(Soil = str_extract(pot, "\\w"),
         ET = potmass2 - potmass3,
         soil_LSB = case_when(Soil == "L" ~ "L",
                          Soil == "S" ~ "S",
                          Soil == "C" ~ "B",
                          Soil == "N" ~ "B"),
         plant_YN = case_when(Soil == "L" ~ "Y", 
                           Soil == "S" ~ "Y",
                           TRUE ~ "N"),
         batch_discharge = case_when(dis_table == 2 ~ 1,
                           dis_table == 3 ~ 1,
                           dis_table == 4 ~ 1,
                           TRUE ~ 2)) %>%
  dplyr::select(pot,ET,Dmass_tot, soil_LSB, plant_YN, batch_discharge, dis_table)

#condense repeat measurements into averages
root1 <- root %>% group_by(pot) %>%
  summarise(mean_rootlength = median(root_cm))

#start with LiCor data because this may have the most measurements
dat1 <- merge(dat, licor, by = "pot", all = T)
dat2 <- merge(dat1, bm, by = "pot", all.x = T)
dat3 <- merge(dat2, pw_geo, by = "pot", all.x = T)
dat4 <- merge(dat3, dis1, by = "pot", all.x = T)
dat5 <- merge(dat4, s_geo, by = "pot", all.x = T)
dat6 <- merge(dat5, root1, by = "pot", all.x = T)

dat7 <- dat6 %>%
  #removes random "X" and "notes" columns
  dplyr::select(1:12, 14:21, 24:55, 57:64) %>%
  #calculate final biomass columns
  mutate(totBM = massA + massB,
         BMrat = massB/massA) 

#Export the finished dataframes ------------------------------------------------
write.csv(census1, file = "data/ModData/FullCensusTimeSeries.csv",row.names = F)
write.csv(dat7, file = "data/ModData/AllPerformanceGeochem.csv", row.names = F)
