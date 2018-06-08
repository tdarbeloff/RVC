#Checking to see if WML is correlated with what it should be, as a sanity check

library(tidyverse)
library(car)
library(ggplot2)
library(haven)
library(jtools)
library(stats)
library(ggpubr)
library(psych)
library(lm.beta)
library(investr)
library(MASS)
library(lmSupport)
library(caret)


d_wml <- read_csv("H:/Database/DBIS/Imaging/WhiteMatterHyperintensities.csv")
d_wml$DWMHvol_mm3[which(is.nan(d_wml$DWMHvol_mm3))]=NA
d_wml <- na.omit(d_wml)

d_rvc38 <- read_sav("Phase38_varsSES.sav")
d_rvc38 <- d_rvc38 %>%
  mutate(gsub("[[:space:]]", "", d_rvc38$subID),
         to_add = ifelse(nchar(d_rvc38$subID) == 3, "DMHDS0", "DMHDS"),
         ID = paste(to_add, subID, sep = ""))
d_rvc38 <- na.omit(d_rvc38)

d_BP <- read_sav("Bloodpressure45_Batch27.sav")
d_BP <- d_BP %>%
  mutate(gsub("[[:space:]]", "", d_BP$subID),
         to_add = ifelse(nchar(d_BP$subID) == 3, "DMHDS0", "DMHDS"),
         ID = paste(to_add, subID, sep = ""))


d_rvc <- as.data.frame(read_sav("Retinal45_combined_wk64.sav"))

d_rvc <- d_rvc %>%
  mutate(gsub("[[:space:]]", "", d_rvc$subID),
         to_add = ifelse(nchar(d_rvc$subID) == 3, "DMHDS0", "DMHDS"),
         ID = paste(to_add, subID, sep = ""))
d_rvc <- na.omit(d_rvc)

#now that datasets area read in, join them together by subID
d <- inner_join(d_BP, d_wml, by=c("ID"))
d <- inner_join(d, d_rvc, by=c("ID"))
d <- inner_join(d, d_rvc38, by=c("ID"))

#z-score retinal data
d$Big6craep45_z <- scale(as.numeric(as.character(d$Big6craep45)))[,1]
d$Big6crvep45_z <- scale(as.numeric(as.character(d$Big6crvep45)))[,1]
d$Big6craep38_z <- scale(as.numeric(as.character(d$Big6craep38)))[,1]
d$Big6crvep38_z <- scale(as.numeric(as.character(d$Big6crvep38)))[,1]

#Transform wmh by log and z-score (retinal vessel caliber has already been checked and found normal)
d$wholeBrainWMHvol_mm3_lg <- log(as.numeric(as.character(d$wholeBrainWMHvol_mm3)))
d$wholeBrainWMHvol_mm3_z <- scale(d$wholeBrainWMHvol_mm3_lg)[,1]

d$PVWMHvol_mm3_z <- scale(log(as.numeric(as.character(d$PVWMHvol_mm3))))[,1]
d$DWMHvol_mm3_z <- scale(log(as.numeric(as.character(d$DWMHvol_mm3))))[,1]

#transform and scale BP variables as necessary
d$MAP45_z <- scale(as.numeric(as.character(d$MAP45)))[,1]
d$Sysbp45_z <- scale(as.numeric(as.character(d$Sysbp45)))[,1]
d$Diabp45_z <- scale(as.numeric(as.character(d$Diabp45)))[,1]


#now lets check correlations between WMH and BP, etc.
cor.test(d$MAP45_z, d$wholeBrainWMHvol_mm3_z)
cor.test(d$MAP45_z, d$Big6craep45_z)
cor.test(d$MAP45_z, d$Big6crvep45_z)
cor.test(d$Sysbp45_z, d$wholeBrainWMHvol_mm3_z)
 

