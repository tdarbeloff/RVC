
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


#Call and format Datasets
d_rvc <- as.data.frame(read_sav("Retinal45_combined_wk64.sav"))

d_rvc <- d_rvc %>%
  mutate(gsub("[[:space:]]", "", d_rvc$subID),
         to_add = ifelse(nchar(d_rvc$subID) == 3, "DMHDS0", "DMHDS"),
         ID = paste(to_add, subID, sep = ""))
d_rvc <- na.omit(d_rvc)

d_rvc38 <- as.data.frame(read_sav("Phase38_vars.sav"))
d_rvc38 <- d_rvc38 %>%
  mutate(gsub("[[:space:]]", "", d_rvc38$subID),
         to_add = ifelse(nchar(d_rvc38$subID) == 3, "DMHDS0", "DMHDS"),
         ID = paste(to_add, subID, sep = ""))
d_rvc38 <- na.omit(d_rvc38)

d_wml <- read_csv("N:/DBIS.01/Data/ALL_DATA_TO_USE/Imaging/WhiteMatterHyperintensities.csv")
d_wml$DWMHvol_mm3[which(is.nan(d_wml$DWMHvol_mm3))]=NA
d_wml <- na.omit(d_wml)

d <- inner_join(d_rvc, d_wml, by=c("ID"))

d <- inner_join(d, d_rvc38, by=c("ID"))
Check Variables
hist(d$Big6crvep45)
summary(d$Big6crvep45)
d$VenCal_mean <- mean(d$Big6crvep45)
d$VenCal_sd <- sd(d$Big6crvep45)
d$VenCal_outlier <- d$Big6crvep45 > (d$VenCal_mean + 3*d$VenCal_sd)
d %>%
  filter(d$VenCal_outlier == TRUE)

#z-score for interpretable/comparable betas
d$Big6crvep45_z <- scale(d$Big6crvep45)[,1]
Big6Venular has normal distribution and no outliers. As of right now, I am not going to transform it in any way.

hist(d$Big6craep45)
summary(d$Big6craep45)
d$ArtCal_mean <- mean(d$Big6craep45)
d$ArtCal_sd <- sd(d$Big6craep45)
d$ArtCal_outlier <- d$Big6craep45 > (d$ArtCal_mean + 3*d$ArtCal_sd)
d %>%
  filter(d$ArtCal_outlier == TRUE)
#scale
d$Big6craep45_z <- scale(d$Big6craep45)[,1]

#Big6ArtCal looks fairly normal (although not quite as good as ven distrib), so as of right now I am not going to transform it, but I may come back to it. There are no outliers in the main distribution.

#Both residual distributions should be normal, but I will check them over anyway.

hist(d$VenularRsd45)
hist(d$ArteriolRsd45)
both look fairly normal. Venular is slightly positively skewed, Arteriole is slightly negatively skewed. This might be something to take note of later in residual analyses, but I am going to be using the non-resid variables for now, so I will put these aside.

#Now to check main WML variables

d$wholeBrainWMHvol_mm3 <- as.numeric(as.character(d$wholeBrainWMHvol_mm3))
hist(d$wholeBrainWMHvol_mm3)
summary(d$wholeBrainWMHvol_mm3)

d$wholeBrainWMHvol_mm3_lg <- log(d$wholeBrainWMHvol_mm3)
hist(d$wholeBrainWMHvol_mm3_lg)

#z-score for readables betas
d$wholeBrainWMHvol_mm3_z <- scale(d$wholeBrainWMHvol_mm3_lg)[,1]


#Isolate WBV frontal as WML are often seen earlier in age in frontal regions as opposed to occipital, etc.
d$Rfrontal_WMHvol_mm3 <- as.numeric(as.character(d$Rfrontal_WMHvol_mm3))
hist(d$Rfrontal_WMHvol_mm3)
summary(d$Rfrontal_WMHvol_mm3)
#Shift over by 1 so it can be log transformed
d$Rfrontal_WMHvol_mm3 <- d$Rfrontal_WMHvol_mm3 + 1
d$Rfrontal_WMHvol_mm3_lg <- log(d$Rfrontal_WMHvol_mm3)
hist(d$Rfrontal_WMHvol_mm3_lg)

d$Lfrontal_WMHvol_mm3 <- as.numeric(as.character(d$Lfrontal_WMHvol_mm3))
hist(d$Lfrontal_WMHvol_mm3)
d$Lfrontal_WMHvol_mm3 <- d$Lfrontal_WMHvol_mm3 + 1
d$Lfrontal_WMHvol_mm3_lg <- log(d$Lfrontal_WMHvol_mm3)
hist(d$Lfrontal_WMHvol_mm3_lg)

#very skewed, but log transforms it into a beautiful normal distrib.

d$wholeBrain_WMHnoc_total <- as.numeric(as.character(d$wholeBrain_WMHnoc_total))
hist(d$wholeBrain_WMHnoc_total)
d$wholeBrain_WMHnoc_total_lg <- log(d$wholeBrain_WMHnoc_total)
hist(d$wholeBrain_WMHnoc_total_lg)
summary(d$wholeBrain_WMHnoc_total_lg)

d$wholeBrain_WMHnoc_total_z <- scale(d$wholeBrain_WMHnoc_total_lg)[,1]


#again, needed the log transformation

#that looks fine, lets move on. Im going to split whole brain WML volume into PVMWH and DWMH, as they have been associated with separate things (PV with more cognitive, etc and D with depression, etc). So lets look at those variables

d$PVWMHvol_mm3 <- as.numeric(d$PVWMHvol_mm3)
hist(d$PVWMHvol_mm3)
d$PVWMHvol_mm3_lg <- log(d$PVWMHvol_mm3)
hist(d$PVWMHvol_mm3_lg)
summary(d$PVWMHvol_mm3)

d$PVWMHvol_mm3_z <- scale(d$PVWMHvol_mm3_lg)[,1]
#log transformed.

d$DWMHvol_mm3 <- as.numeric(d$DWMHvol_mm3)
hist(d$DWMHvol_mm3)
summary(d$DWMHvol_mm3)
#add 1 to shift over so we can log transform
d$DWMHvol_mm3 <- d$DWMHvol_mm3 + 1
summary(d$DWMHvol_mm3)
d$DWMHvol_mm3_lg <- log(d$DWMHvol_mm3)
hist(d$DWMHvol_mm3_lg)

d$DWMHvol_mm3_z <- scale(d$DWMHvol_mm3_lg)[,1]

#That log transformation left a slight negative skew. I am going to leave it for right now, but I may come back and work on this one more.

#I am going to take a peek at one of the tortuosity variables, however I need to know more before I wade into all of them.

hist(d$cTORTv_p45)
summary(d$cTORTv_p45)
d$cTORTv_p45_lg <- log(d$cTORTv_p45)
hist(d$cTORTv_p45_lg)

d$cTORTv_p45_z <- scale(d$cTORTv_p45_lg)[,1]
#I also want to look at IQ (age 38), so lets look at those.

hist(d$WFSIQ38)
summary(d$WFSIQ38)
d$WFSIQ38_z <- scale(d$WFSIQ38)

##Now to try some preliminary regressions

#we will start with ven and art (major 6) regressed on whole brain WML volume

wb_ven <- lm(wholeBrainWMHvol_mm3_z ~ Big6crvep45_z, d)
summary(wb_ven)

wb_art <- lm(wholeBrainWMHvol_mm3_z ~ Big6craep45_z, d)
summary(wb_art)

#SO from preliminary analyses (no exclusions in data, etc), retinal venular caliber is not associated with WML, but arteriole is (very small effect).

#lets split whole brain by PVWML

pv_ven <- lm(PVWMHvol_mm3_z ~ Big6crvep45_z, d)
summary(pv_ven)

pv_art <- lm(PVWMHvol_mm3_z ~ Big6craep45_z, d)
summary(pv_art)


#similar result as above. Lets check if the trend holds with DWMH. Im guessing that neither will be significant

d_ven <- lm(DWMHvol_mm3_z ~ Big6crvep45_z, d)
summary(d_ven)

d_art <- lm(DWMHvol_mm3_z ~ Big6craep45_z, d)
summary(d_art)


#as expected, neither are significantly associated.

#lets look to replicate PVWMH associated with decreased cognition

pv_cog <- lm(PVWMHvol_mm3_z ~ WFSIQ38_z, d)
summary(pv_cog)

#that replicates, which is nice. And in a younger cohort than probably ever seen. Which is also nice. This is definitely something we should follow up on.

#What happens if we move to an earlier timepoint? (RVC age 38)

#look at variables quick
hist(d$Big6craep38)
hist(d$Big6crvep38)

d$Big6craep38_z <- scale(d$Big6craep38)[,1]
d$Big6crvep38_z <- scale(d$Big6crvep38)[,1]

summary(lm(wholeBrainWMHvol_mm3_z ~ Big6craep38_z, d))
summary(lm(wholeBrainWMHvol_mm3_z ~ Big6crvep38_z, d))
#so RVC at age 38 does not predict WML volume.

#brief residual analysis (one of many to come) to double check linear regression is the right move here.

d$resid_pvArt <- pv_art$residuals
d$predicted_pvArt <- pv_art$fitted.values
hist(d$resid_pvArt)

ggplot(d, aes(y=resid_pvArt, x=predicted_pvArt, label = rownames(d))) + geom_point() +
geom_smooth(method="loess") + geom_smooth(method="lm", se=FALSE, color="red") +
geom_text(hjust=0, nudge_x=.05)
avPlots(pv_art)