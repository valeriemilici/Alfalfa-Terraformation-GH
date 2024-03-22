### This sheet will briefly summarize the biomass of the potting soil controls

rm(list =ls())

library(tidyverse)

dat <- read.csv("data/ModData/pottingsoilmass.csv")


## Split data so that type (A or B) is two columns for aboveground and
## belowgrownd biomass

dat1 <- pivot_wider(data = dat, names_from = "type", values_from = "mass") %>%
  rename(MassA = A, MassB = B) %>%
  mutate(MassT = MassA + MassB)

#What's the range of masses?

range(dat1$MassT)
#0.23 - 0.94
quantile(dat1$MassT)
#0.47

dat2 <- read.csv("data/ModData/AllPerformanceGeochem.csv")

t.test(dat1$MassT, dat2$totBM)
