### An analysis to generate some PCAs to see if the genotypes cluster around
### various traits and response variables.

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(corrr) #for PCAs
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggfortify)
library(ggplot2)

dat <- read.csv("data/performance.csv") 
datdis <- read.csv("data/discharge.csv") #get et for each pot
alfdat <- read.csv("data/GenoSelections/AlfPISelections.csv") #alfalfa data
geochem <- read.csv("data/disgeochem.csv")
#### Data Cleaning -------------------------------------------------------------

#prep alfalfa trait data for merging with master dataset
Geno <- c("D16", "A20", "G15", "212", "SWE" ,"R24", "TAS", "YON", "K19",
          "MOC", "INA", "CRI", "639", "105", "G74", "515", "VIR")
geno <- c("other", "other", "G15", "other", "other", "other", "TAS", "YON", 
          "other", "MOC", "other", "CRI", "other", "other", "other", "other", "VIR")

alfdat2 <- cbind(Geno, alfdat)

alfdat3<- alfdat2[,c(1,6,7,10,11)]

#keep only the columns that we need
datdis1 <- datdis %>% dplyr::select(pot, et)
dat0 <- left_join(dat, datdis1, by = "pot")
dat00 <- left_join(dat0, alfdat3, by = "Geno")
dat1 <- dat00 %>% filter(CensusNo == 6,
                       !is.na(RGR_d), !is.na(meanTPN),
                       !is.na(massA), !is.na(et)) %>%
  dplyr::select(Height, LeafNo, massA, massB,rootsample, meanPSN, meanTPN,
         Soil, Geno, RGR_d, et, Days.Plant.to.Maturity,
         Determinate.Taproot.Percentage, Fibrous.Root.Mass, pot)

colSums(is.na(dat1))
#good
pot <- dat1$pot

#make sure your data is only numbers
dat2 <- dat1 %>% dplyr::select(!c(Geno, rootsample, Soil, pot))

### https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
Geno <- data.frame(Geno)
geno<- data.frame(geno)
type <- cbind(Geno, geno)

dat1b<- left_join(dat1, type, by = "Geno")

#create the PCA and have the function internally scale the variables
perf_pca <- PCA(dat2)
#let's make sure we have a good description of each PC

eig.val <- get_eigenvalue(pca_res)


#we have 11 variables in the PCA, 100/11 = 9.09, so any eigenvalue that is 
#greater than 9.09 indicates a PC that explains more variance than a single
#variable alone. 
View(eig.val)

#Dims 1-3 all explain more than 9.09 of the variance. By a non-arbitrary 
#cutoff we should use 1-3. I will move forward with 1 & 2 only just to teach
#myself, but will come back to dim 3 at a later point once I know more. 

#What components dominate PC1?
fviz_contrib(pca_res, choice = "var", axes = 1, top = 10)
#Height, massA&B, LeafNo, and MeanPSN. All traits associated with growth.

fviz_contrib(pca_res, choice = "var", axes = 2, top = 10)
#Transpiration rate and days to maturity (a little bit PSN). This is a bit more
#challenging to define. Water loss and growth rates...

#because PC3 did come out as important, do the work now
fviz_contrib(pca_res, choice = "var", axes = 3, top = 10)
#fibrous root mass, et, PSN, and days to maturity. Stuff relating to how 
#root morphology affects water usage and is maybe affected by how long it takes
#a plant to reach maturity? (note- none of our plants actually reached a mature
#growth form. 

res.desc <- dimdesc(perf_pca, axes = c(1,2), proba = 0.05)
res.desc$Dim.1 #all vars are positively correlated
res.desc$Dim.2 #days to maturity is pos cor, but TPN is negative. So high TPN
#and fast growth are related, or they trade off?

#plot the PCA, the "geno" column has all non-interesting genos listed as "other"
autoplot(pca_res,
         data = dat1b,
         color = "geno",
         loadings = T,
         loadings.label = T) +
  theme_bw()

#save all of the individual point data (coords and loadings) to a new dataframe
ind <- get_pca_ind(pca_res)

#the coordinates on all PCs for each data point. This information can be used
#in follow-up models to assess if a constellation of performance data 
#correlates to any genos
perf_coords_all <- cbind(data.frame(ind$coord[,1:2]),pot)
write.csv(perf_coords_all, "data/perfcoords2.csv")

### Performance PCA no root traits ---------------------------------------------

#[1] get the data ready ----------

dat3 <- dat2 %>% select(!c(Days.Plant.to.Maturity, Determinate.Taproot.Percentage, 
                           Fibrous.Root.Mass))

#[2] Make PCA ------------------
perf_pca_2 <- PCA(dat3)

#[3] Describe dimentions --------

get_eigenvalue(perf_pca_2)
#Dims 1-3 once more could be useful, but we'll focus on 1 & 2

fviz_contrib(perf_pca_2, choice = "var", axes = 1)
#defined as size-based traits (height, mass, leaf number)

fviz_contrib(perf_pca_2, choice = "var", axes = 2)
#gas-exchange traits: PSN and TPN

fviz_contrib(perf_pca_2, choice = "var", axes = 3)
#pure growth rate...

PCA.desc <- dimdesc(perf_pca_2, axes = 1:2, proba = 0.05)
PCA.desc$Dim.1
#all things are positively correlated

PCA.desc$Dim.2
#TPN and PSN are positive correlated 

#[4] plot it ------------------
fviz_pca_biplot(perf_pca_2, 
                col.ind = dat1b$geno, palette = "jco", 
                addEllipses = T, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Geno") 

#[5] Extract the data that you need --------

ind2 <- get_pca_ind(perf_pca_2)

perf_coords_response <- cbind(data.frame(ind2$coord[,1:2]), pot)
write.csv(perf_coords_response, "data/perfcoords1.csv")

### Discharge Geochemistry PCA -------------------------------------------------

#[1] get the data ready ----------
chem2 <- geochem[,c(2,4,8:21)] 
ins <- c("MOC", "YON", "CRI", "TAS", "VIR")
chem2 <- chem2 %>% filter(Geno %in% ins)

chem2 <- chem2 %>% mutate(geno = case_when(Geno == "MOC" ~ "MOC",
                                           Geno == "CRI" ~ "CRI",
                                           Geno == "YON" ~ "YON",
                                           Geno == "TAS" ~ "TAS",
                                           Geno == "VIR" ~ "VIR",
                                           TRUE ~ "other"))

colSums(is.na(chem2)) #remove nitrite, phosphate, and lithium bc so many NAs
#remove Carbon bc we don't have control data yet, and would remove all control
chem3 <- chem2 %>% dplyr::select(!c(Nitrite, Phosphate, Lithium, Carbon, Nitrate,
                             Calcium)) %>%
  filter(!is.na(Fluorine), !is.na(Chlorine), !is.na(Sulfate),
         !is.na(pH), !is.na(EC), !is.na(Sodium), !is.na(Potassium),
         !is.na(Magnesium) #, !is.na(Calcium)
         )

colSums(is.na(chem3))

chem4 <- chem3[,-c(1,2,11)]
#[2] Make PCA ------------------
dis_pca <- PCA(chem4)

#[3] Describe dimensions --------

get_eigenvalue(dis_pca)

fviz_contrib(dis_pca, choice = "var", axes = 1)
#Sodium, Potassium,Magnesium

fviz_contrib(dis_pca, choice = "var", axes = 2)
#Fluorine, Chlorine, Sulfate

PCA.desc <- dimdesc(dis_pca, axes = 1:2, proba = 0.05)
PCA.desc$Dim.1
#all things are positively correlated

PCA.desc$Dim.2
#all main elements are postively correlated

#[4] plot it ------------------
porewater <- fviz_pca_ind(dis_pca, 
                col.ind = chem3$Geno, palette = "jco", 
                addEllipses = T, ellipse.level = 0.9, label = "none", repel = TRUE,
                legend.title = "Geno") +
  theme_classic(12)

porewater

ggsave(plot = porewater, filename = "figures/porewaterPCA.png")
#[5] Extract the data that you need --------

ind3 <- get_pca_ind(dis_pca)

dis_coords <- cbind(data.frame(ind3$coord[,1:2]), chem3$pot)
write.csv(dis_coords, "data/discoords.csv")

### Geochem and Performance PCA ------------------------------------------------

#[1] get the data ready ----------

dat4 <- anion %>% group_by(pot) %>%
  summarise(Fluorine = median(F.mg.L.),
            Nitrate = median(NO3.mg.L.),
            Sulfate = median(SO4.mg.L.),
            Chlorine = median(Cl.mg.L.)) %>%
  filter(!is.na(Fluorine),
         !is.na(Nitrate),
         !is.na(Sulfate),
         !is.na(Chlorine)) %>%
  mutate(Geno = str_extract(pot, "\\w{2,}"))

dat5 <- left_join(dat4, pH, by = "pot")
dat5 <- left_join(dat5, type, by = "Geno")

dat7<- left_join(dat5, dat1, by = "pot")

dat7 <- filter(dat7, !is.na(Height))
dat8 <- dat7 %>% select(!c(flat, pot, X, volume..ml., geno, Fibrous.Root.Mass,
                           Determinate.Taproot.Percentage, rootsample, Geno.x,
                           Soil, Geno.y, Days.Plant.to.Maturity))

dat8$Chlorine <- as.numeric(dat8$Chlorine)

#[2] Make PCA ------------------
all_pca <- PCA(dat8)

#[3] Describe dimentions --------

get_eigenvalue(all_pca)
#First 5 PCs all have some power

fviz_contrib(all_pca, choice = "var", axes = 1)
#Chlorine, Sulfate, et, pH

fviz_contrib(all_pca, choice = "var", axes = 2)
#MassA&B, Height, LeafNO, PSN (AKA size stuff), and PH

PCA.desc <- dimdesc(all_pca, axes = 1:2, proba = 0.05)
PCA.desc$Dim.1
#Geochem is all positive, et is negative

PCA.desc$Dim.2
#All major contributors are positively correlated

#[4] plot it ------------------
fviz_pca_biplot(all_pca, 
                col.ind = dat7$geno, palette = "jco", 
                addEllipses = F, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Geno") 

#[5] Extract the data that you need --------

ind4 <- get_pca_ind(all_pca)

all_coords <- cbind(data.frame(ind4$coord[,1:2]), dat7$pot)
write.csv(all_coords, "data/allcoords.csv")
#this one might have too little information (only 43 pots had all info)

###
