### Bootstrap the model output

### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(lme4) #mixed effects models
library(lmerTest) #p-values from linear mixed models
library(ggplot2) # plot results
library(patchwork) #joins plots nicely
library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

### Read in Models -------------------------------------------------------------
#plant performance
TotBM <- readRDS("04_analyses/02_models/output/BMSoilT.RDS")
RootPropBM <- readRDS("04_analyses/02_models/output/BMSoilR.RDS")
Surv <- readRDS("04_analyses/02_models/output/SurvGeno1.RDS")
#Geochemistry Soil
EC <- readRDS("04_analyses/02_models/Output/ECSoil.RDS")
pH <- readRDS("04_analyses/02_models/Output/pHSoil.RDS")
Carbon <- readRDS("04_analyses/02_models/Output/CSoil.RDS")
Calcium <- readRDS("04_analyses/02_models/Output/CaSoil.RDS")
Potassium <- readRDS("04_analyses/02_models/Output/KSoil.RDS")
Magnesium <- readRDS("04_analyses/02_models/Output/MgSoil.RDS")
Sodium <- readRDS("04_analyses/02_models/Output/NaSoil.RDS")
Nitrate <- readRDS("04_analyses/02_models/Output/NO3Soil.RDS")
IC <- readRDS("04_analyses/02_models/Output/ICSoil.RDS")
#Geochemistry Geno
ECGeno <- readRDS("04_analyses/02_models/Output/ECGeno.RDS")
CGeno <- readRDS("04_analyses/02_models/Output/CGeno.RDS")
CaGeno <- readRDS("04_analyses/02_models/Output/CaGeno.RDS")
KGeno <- readRDS("04_analyses/02_models/Output/KGeno.RDS")
MgGeno <- readRDS("04_analyses/02_models/Output/MgGeno.RDS")
NaGeno <- readRDS("04_analyses/02_models/Output/NaGeno.RDS")
NO3Geno <- readRDS("04_analyses/02_models/Output/NO3Geno.RDS")
ICGeno <- readRDS("04_analyses/02_models/Output/ICGeno.RDS")

# I will now create a series of functions where I simulate data from the model
# output and then bootstrap the model predictions to extract the CIs. I will
# number each of these as 1-7 according to the order in which the models are 
# listed above.

### Boostrap Function ----------------------------------------------------------
##Performance ----------------------------------------
#[1]
BR.fun <- function(.) {
  preddat <- expand.grid(Soil = "L",
                         Level =c("1", "2", "3"),
                         totBM = 0.067,
                         Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))
  predict(., newdata = preddat, re.form = ~0)
}

#[2]
BT.fun <- function(.) {
  preddat <- expand.grid(SeedRange = c("Native", "Invaded"),
                      SoilRangeInvaded = c("In", "FB", "Out"),
                      Repromass = 21)
  predict(., newdata = preddat, re.form = ~0)
}

#[3]
BST.fun <- function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         Level =c("1", "2", "3"))
  predict(., newdata = preddat, re.form = ~0)
}

#[4]
BSR.fun <- function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         totBM = 0.089,
                         Level =as.factor(2))
  predict(., newdata = preddat, re.form = ~0)
}

#[5]
PSN.fun <- function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         Level =as.factor(2),
                         Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))
  predict(., newdata = preddat, re.form = ~0)
}

#[6]
TPN.fun <- function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         Level =as.factor(2),
                         Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))
  predict(., newdata = preddat, re.form = ~0)
}

#[7]
Surv.fun <- function(.) {
  preddat <- expand.grid(Soil = c("L", "S"),
                         Level =as.factor(2),
                         Geno = c("MOC", "TAS", "CRI", "VIR", "YON"))
  predict(., newdata = preddat, re.form = ~0)
}

##Geochemistry ------------------------------------------
#[8]
geochem.fun <- function(.) {
  preddat <- expand.grid(Soil = c("C", "N", "S", "L"),
                         Level =c("1", "2", "3"),
                         Dmass_tot = 106)
  predict(., newdata = preddat, re.form = ~0)
}

#[9]
C.fun <- function(.) {
  preddat <- expand.grid(Soil = c("C", "N", "S", "L"),
                         Level =c("1", "2", "3"),
                         Dmass_tot = 106,
                         TC_batch = "1")
  predict(., newdata = preddat, re.form = ~0)
}

#[10] 
geochem.geno.fun <- function(.) {
preddat <- expand.grid(Soil = c("S", "L"),
                       Level =c("1", "2", "3"),
                       Geno = c("CRI", "K19", "MOC",
                                "TAS", "VIR", "G15", "YON"),
                       Dmass_tot = 106)
predict(., newdata = preddat, re.form = ~0)
}

#[11]
cation.fun <- function(.) {
  preddat <- expand.grid(Soil = c("S", "L"),
                         Level =c("1", "2", "3"),
                         Dmass_tot = 106,
                         Geno = c("MOC", "VIR", "G15", "105"))
  predict(., newdata = preddat, re.form = ~0)
}

#[12]
IC.fun <- function(.) { #note this function isn't working for an unknown reason
  preddat <- expand.grid(Soil = c("C", "N", "S", "L"),
                         Level =c("1", "2", "3"),
                         Dmass_tot = 106,
                         batch = c("A", "B", "C"))
  predict(., newdata = preddat, re.form = ~0)
}

### Bootstrap the Model --------------------------------------------------------
##Performance -------------------------------------------
#[1]
BR.boot <- bootMer(BioRat, nsim = 1000, FUN = BR.fun,
                         parallel="snow", ncpus = detectCores(), 
                         cl=cl)

#[2]
BT.boot <- bootMer(BioTot, nsim = 1000, FUN = BT.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)

#[3]
BST.boot <- bootMer(BmSoilTot, nsim = 1000, FUN = BST.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)

#[4]
BSR.boot <- bootMer(BmSoilRat, nsim = 1000, FUN = BSR.fun,
                    parallel="snow", ncpus = detectCores(), 
                    cl=cl)

#[5]
PSN.boot <- bootMer(PSN, nsim = 1000, FUN = PSN.fun,
                    parallel="snow", ncpus = detectCores(), 
                    cl=cl)

#[6]
TPN.boot <- bootMer(TPN, nsim = 1000, FUN = TPN.fun,
                    parallel="snow", ncpus = detectCores(), 
                    cl=cl)

#[7]
Surv.boot <- bootMer(Surv, nsim = 1000, FUN = Surv.fun,
                    parallel="snow", ncpus = detectCores(), 
                    cl=cl)
##Geochemistry Soil --------------------------------------
#[8]
EC.boot <- bootMer(EC, nsim = 1000, FUN = geochem.fun,
                     parallel="snow", ncpus = detectCores(), 
                     cl=cl)
#[9]
pH.boot <- bootMer(pH, nsim = 1000, FUN = geochem.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[10]
C.boot <- bootMer(Carbon, nsim = 1000, FUN = C.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[11]
Ca.boot <- bootMer(Calcium, nsim = 1000, FUN = cation.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[12]
Mg.boot <- bootMer(Magnesium, nsim = 1000, FUN = cation.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[13]
NO3.boot <- bootMer(Nitrate, nsim = 1000, FUN = geochem.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[14]
K.boot <- bootMer(Potassium, nsim = 1000, FUN = cation.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[15]
Na.boot <- bootMer(Sodium, nsim = 1000, FUN = cation.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[15B]
IC.boot <- bootMer(IC, nsim = 1000, FUN = IC.fun,
                  parallel="snow", ncpus = detectCores(), 
                  cl=cl)

## Geochemistry Geno --------------------------------------
#[16]
ECGeno.boot <- bootMer(ECGeno, nsim = 1000, FUN = geochem.geno.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)
#[17]
CGeno.boot <- bootMer(CGeno, nsim = 1000, FUN = geochem.geno.fun,
                       parallel="snow", ncpus = detectCores(), 
                       cl=cl)
#[18]
CaGeno.boot <- bootMer(CaGeno, nsim = 1000, FUN = geochem.geno.fun,
                       parallel="snow", ncpus = detectCores(), 
                       cl=cl)
#[19]
KGeno.boot <- bootMer(KGeno, nsim = 1000, FUN = geochem.geno.fun,
                       parallel="snow", ncpus = detectCores(), 
                       cl=cl)
#[20]
MgGeno.boot <- bootMer(MgGeno, nsim = 1000, FUN = geochem.geno.fun,
                       parallel="snow", ncpus = detectCores(), 
                       cl=cl)
#[21]
NaGeno.boot <- bootMer(NaGeno, nsim = 1000, FUN = geochem.geno.fun,
                       parallel="snow", ncpus = detectCores(), 
                       cl=cl)
#[22]
NO3Geno.boot <- bootMer(NO3Geno, nsim = 1000, FUN = geochem.geno.fun,
                       parallel="snow", ncpus = detectCores(), 
                       cl=cl)

#[23]
ICGeno.boot <- bootMer(ICGeno, nsim = 1000, FUN = geochem.geno.fun,
                      parallel="snow", ncpus = detectCores(), 
                      cl=cl)
#stop the clusters once bootstrapping is finished-------
stopCluster(cl = cl) 

### Save the Output ------------------------------------------------------------
##Performance ---------------------------------------
save(BR.boot, file = "04_Analyses/03_Bootstrap/output/BR.boot")
save(BT.boot, file = "04_Analyses/03_Bootstrap/output/BT.boot")
save(BSR.boot, file = "04_Analyses/03_Bootstrap/output/BSR.boot")
save(BST.boot, file = "04_Analyses/03_Bootstrap/output/BST.boot")
save(PSN.boot, file = "04_Analyses/03_Bootstrap/output/PSN.boot")
save(TPN.boot, file = "04_Analyses/03_Bootstrap/output/TPN.boot")
save(Surv.boot, file = "04_Analyses/03_Bootstrap/output/Surv.boot")
##Geochemistry Soil PW ---------------------------------------
save(EC.boot, file = "04_Analyses/03_Bootstrap/output/EC.boot")
save(pH.boot, file = "04_Analyses/03_Bootstrap/output/pH.boot")
save(C.boot, file = "04_Analyses/03_Bootstrap/output/C.boot")
save(Ca.boot, file = "04_Analyses/03_Bootstrap/output/Ca.boot")
save(Mg.boot, file = "04_Analyses/03_Bootstrap/output/Mg.boot")
save(NO3.boot, file = "04_Analyses/03_Bootstrap/output/NO3.boot")
save(K.boot, file = "04_Analyses/03_Bootstrap/output/K.boot")
save(Na.boot, file = "04_Analyses/03_Bootstrap/output/Na.boot")
save(IC.boot, file = "04_analyses/03_Bootstrap/output/IC.boot")
##Geochemistry Geno --------------------------------------------
save(ECGeno.boot, file = "04_Analyses/03_Bootstrap/output/ECGeno.boot")
save(CGeno.boot, file = "04_Analyses/03_Bootstrap/output/CGeno.boot")
save(CaGeno.boot, file = "04_Analyses/03_Bootstrap/output/CaGeno.boot")
save(MgGeno.boot, file = "04_Analyses/03_Bootstrap/output/MgGeno.boot")
save(KGeno.boot, file = "04_Analyses/03_Bootstrap/output/KGeno.boot")
save(NaGeno.boot, file = "04_Analyses/03_Bootstrap/output/NaGeno.boot")
save(NO3Geno.boot, file = "04_Analyses/03_Bootstrap/output/NO3Geno.boot")
save(ICGeno.boot, file = "04_Analyses/03_Bootstrap/output/ICGeno.boot")
