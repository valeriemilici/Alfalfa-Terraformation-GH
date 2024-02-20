### Super PCA for FRES proposal
#Includes both geochemistry and performance information

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(corrr) #for PCAs
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(ggplot2)

datdis <- read.csv("data/discharge.csv") #performance data 
geochem <- read.csv("data/disgeochem.csv") #pw geochem data
catdat <- read.csv("data/porewater_cations_2024.csv") #updated pw cation data
NPOC <- read.csv("data/NPOC.csv") #organic carbon- more full of a read?
### Data prep
#remove bare soil controls & unnecessary columns from dataset
catdat <- catdat %>% drop_na(Geno)
datdis2 <- datdis %>% drop_na(Geno) %>%
  select(c(4, 22, 29, 36:38, 40, 44)) 
geochem2 <- geochem %>% drop_na(Geno) %>%
  select(c(2,4,8:15, 21))

#retain Cation readings for MOC, G15, and VIR
ins <- c("MOC", "G15", "VIR")
geochem3 <- geochem %>% drop_na(Geno) %>% 
  mutate(lithium = Lithium,
         sodium = Sodium,
         potassium = Potassium,
         magnesium = Magnesium,
         calcium = Calcium) %>%
  filter(Geno %in% ins) %>%
  select(c(2,4,8:15, 21:26))

addcat <- left_join(geochem3, datdis2, by = "pot")
addcat <- addcat %>% rename(Geno = Geno.x)

#join data together
perf_geo <- left_join(datdis2, geochem2, by = "pot")
perf_geo <- left_join(perf_geo, catdat, by = "pot")

perf_geo2 <- perf_geo %>% filter(!Geno.x %in% ins) #remove faulty genos

#create "geno" column that allows us to focus on Genos of interest
Genos <- data.frame(cbind( Geno = c("D16", "A20", "G15", "212", "SWE" ,"R24", "TAS", "YON", "K19",
          "MOC", "INA", "CRI", "639", "105", "G74", "515", "VIR"),
          geno = c("other", "other", "G15", "other", "other", "other", "TAS", "YON", 
          "other", "MOC", "other", "CRI", "other", "other", "other", "other", "VIR")))

perf_geo2 <- left_join(perf_geo2, Genos, by = "Geno") #add new column to data
addcat <- left_join(addcat, Genos, by = "Geno")

#rearrange columns in for rbind
addcat2 <- addcat %>% select(1,18:21,23,3:16,2,24)
perf_geo3 <- perf_geo2 %>% select(1,3:6,8,10:18,20:21,23:25, 27,28)

perf_geo4<- rbind(addcat2, perf_geo3)

#add in NPOC column (try to have fuller read of organic carbon)
NPOC <- select(NPOC, 1:2)

perf_geo4 <- left_join(perf_geo4, NPOC, by = "pot")
#filter out the NAs
perf_geo1 <- perf_geo4 %>% filter( !is.na(PSN),
                                !is.na(TPN),
                                !is.na(pH),
                                !is.na(Geno),
                                !is.na(Chlorine),
                                !is.na(NPOC),
                                !is.na(Fluorine),
                                !is.na(sodium))%>%
  select(-c(9:12,15,16,18:20))

colSums(is.na(perf_geo1)) #this required a lot of filtering. I don't love it.

#make sure your data is only numbers
dat <- perf_geo1 %>% select(!c(geno, Geno,pot))

## Create the PCA

#[2] Make PCA ------------------
geoperf_pca <- PCA(dat)

#[3] Describe dimentions --------

get_eigenvalue(geoperf_pca)
#Dims 1-3 once more could be useful, but we'll focus on 1 & 2

fviz_contrib(geoperf_pca, choice = "var", axes = 1)
#defined as some geochem measures: sodium, EC, and pH and then et

fviz_contrib(geoperf_pca, choice = "var", axes = 2)
#height, chlorine, and fluorine

fviz_contrib(geoperf_pca, choice = "var", axes = 3)
#TPN, PSN, height (performance)

PCA.desc <- dimdesc(geoperf_pca, axes = 1:2, proba = 0.05)
PCA.desc$Dim.1
#all things are positively correlated

PCA.desc$Dim.2
#TPN and PSN are positive correlated 

#[4] plot it ------------------
fviz_pca_biplot(geoperf_pca, 
                col.ind = perf_geo1$geno, palette = "jco", 
                addEllipses = T, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Geno") 


         