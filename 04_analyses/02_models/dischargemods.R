### Initialize workspace--------------------------------------------------------
rm(list =ls())

library(tidyverse) #clean data manipulation
library(ggplot2) # plot results
library(lme4)
library(lmerTest)
library(performance)

#read in all performance data
dat <- read.csv("data/discharge.csv")

# Describe discharge per treatment
dat2 <- dat %>% filter(quality != "N"&
                       Dmass_tot <=200) #remove poor-quality samples

dat2$Level <- as.character(dat2$Level)

dis_summary <- ggplot(dat2, aes(Level,Dmass_tot)) +
  geom_boxplot() +
  ylab("Total Discharge (ml)") +
  xlab("Percent of Total Water Holding Capacity") +
  scale_x_discrete(labels = c("1" = "30-50", "2" = "55-75", "3" = "80-100")) +
  theme_classic() +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18))

dis_summary

ggsave(dis_summary, file = "figures/dis_summary.jpeg")
# models -----------------------------------------------------------------------

#remove outliers
dat1 <- dat %>% filter(et <= 25,
                       et >= 5)

#model fit isn't quite right. Consider other distributions. So far you've 
#tried lognormal and gamma and normal. The problem seems to be that the model
# is underpredicted the high values. 

mod1b <- lmer(log(et) ~ soil*as.character(Level) + as.character(batch)
              + (1|dis_table),
              data = dat1)
mod1a <- lmer(log(et) ~ soil + as.character(Level) + as.character(batch)
              + (1|dis_table),
              data = dat1)

check_model(mod1)

summary(mod1a)

anova(mod1b, mod1a) #the interaction isn't useful dAIC = ~3
#use mod1a


library(parallel) #allows bootstrapping in parallel
library(stats) #for quantile extraction

#Create cluster within computer to bootstrap in parallel
cl <- makeCluster(detectCores()) 
clusterEvalQ(cl, library(lme4))

dis.fun <- function(.) {
  preddat <- expand.grid(soil = c("L", "S", "B"),
                         Level =c(1:3),
                         batch = c(1,2))
  predict(., newdata = preddat, re.form = ~0)
}

dis.boot <- bootMer(mod1a, nsim = 1000, FUN = dis.fun,
                   parallel="snow", ncpus = detectCores(), 
                   cl=cl)

#stop the clusters once bootstrapping is finished
stopCluster(cl = cl) 

preddis <- expand.grid(soil = c("L", "S", "B"),
                      Level = c(1:3),
                      batch = c(1,2))

preddis$preds <- predict(mod1a, newdata = preddis, re.form = ~0)

preddis<- data.frame(preddis, confint(dis.boot))

names(preddis)[5:6] <- c("lwr", "upr")

# Create Plot

displot <- ggplot() +
  geom_pointrange(data = filter(preddis, batch == 1),
                  mapping = aes(soil, exp(preds),
                                ymin = exp(lwr), ymax = exp(upr),
                                col = as.factor(Level)),
                  size = 1.2,
                  position = position_dodge(width = 0.4)) +
  labs(x = "Live, Sterile, or Bare Soil", y = "Evapotranspiration (ml/day)",
       color = "Water Holding \nCapacity") +
  scale_fill_viridis_d(option = "cividis", direction = -1,
                       labels = c("30-50%", "55-75%", "80-100%")) +
  scale_color_viridis_d(option = "cividis", direction = -1,
                        labels = c("30-50%", "55-75%", "80-100%")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text  = element_text(size = 18)) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))

displot

      ggsave(displot, filename = "figures/discharge.jpeg",
       height = 6, width =8, units = "in")

mod2 <- glmer(log(et) ~ soil + as.character(Level) + as.character(batch)
              + (1|dis_table),
             data = dat1,
             family = Gamma)
check_model(mod2)
summary(mod2)

