## Luke Ozsanlav-Harris

## Updates: 19/08/2022

## Run Cox proportional hazard models on incubation data to model variaiton in nest failure rates
## Use climatic variables as the explanatories in the model and exmaine effect of sub-population
## Need to load the function from the script CoxPH-additional-functions
 

## packages required
pacman::p_load(tidyverse, data.table, survival, finalfit, 
               survminer, coxme, MuMIn, forestplot, dplyr, ltm)

## read in helper functions needed to format data for counting process coxph models
source("Code/Helper functions/CoxPH-additional-functions.R")




#-----------------------------------------#
#### 1. Read in Averaged Env data sets ####
#-----------------------------------------#

## read in climatic data
Window <- fread("Outputs/Average_Env_data_from_arrival_windows_equalfixes.csv")
Roll <- fread("Outputs/Rolling_Env_data_per_day_breeding_season_evenfixes.csv")


## Add extra columns to Roll data set so will bind later on to incubation data
Roll <- separate(Roll, col = tag_date, into = c("ID", "date"), sep = "_")
Roll$Tag_year <- paste0(Roll$ID, "_", year(as.Date(Roll$date)))
Roll$yday <- yday(as.Date(Roll$date))
Roll$ID <- NULL




#---------------------------------------------#
#### 2. Read in Incubation data and format ####
#---------------------------------------------#

## Read in the incubation data set
Inc_ext <- fread("Outputs/Incubation_attempts_extra_update.csv")
Inc_breeders <- Inc_ext

## add tag_date column with tag ID and incubation start date
Inc_breeders$tag_date <- paste0(Inc_breeders$ID, "_", Inc_breeders$attempt_end)




#-------------------------------------#
#### 3. PCA of climatic variables ####
#-------------------------------------#

## add 10 day env windows from greenland arrival to main data set
Inc_br <- inner_join(Inc_breeders, Window, by = "Tag_year")

## examine correlation between climatic variables first
sub1 <- subset(Inc_br, select = c("Gr10precip", "Gr10NDVI", "Gr10Sn", "Gr10temp", "Gr20precip", "Gr20NDVI", "Gr20Sn", "Gr20temp"))
cor_matrix <- cor(sub1)
## Preip doesnt correlate with any other climatic variable, absolute value of r < 0.15
## NDVI and Snow cover have high correlation, absolute  r > 0.75
## Temp has mid correlation with NDVI and snow cover, absolute r c0.5

## Select the varibale for PCA and scale them
sub_pca <- subset(Inc_br, select = c("Gr10NDVI", "Gr10Sn", "Gr10temp"))
sub_pca <- as.data.frame(scale(sub_pca))

## Carry out the PCA
Clim_PCA <- princomp(sub_pca)
summary(Clim_PCA) # summary of componenets
Clim_PCA$loadings # PCA Loadings

## Plot outputs of the PCA
plot(Clim_PCA)
loadings <- as.data.frame(Clim_PCA$loadings[,1:3])
loadings$Symbol <- row.names(loadings)
loadings <- tidyr::gather(loadings, "Component", "Weight", -Symbol)
ggplot(loadings, aes(x=Symbol, y = Weight)) + 
  geom_bar(stat='identity') + 
  facet_grid(Component ~ ., scales = "free_y")

## plot just the component 1 loadings
loadings %>% 
filter(Component == "Comp.1") %>% 
ggplot(aes(x=Symbol, y = Weight)) + 
  geom_bar(stat='identity') +
  theme_bw() +
  xlab("Climatic Variable") +
  scale_x_discrete(labels= c("NDVI", "Snow cover", "Temperature")) +
  theme(panel.grid.major.x = element_blank(), 
        axis.text=element_text(size=12), axis.title=element_text(size=15, face = "bold"), 
        panel.grid.minor.x = element_blank(), strip.text.x = element_text(size = 13))

## read out this plot
ggsave("Paper Plots/Supp Fig 1- PCA axis 1 loadings.png", 
              width = 24, height = 24, units = "cm")


## First two components explain 95% of the varation, will try and add both of them to binary models
## Add first and second principal component to data set
Inc_br$Comp1 <- Clim_PCA$scores[,1]
Inc_br$Comp2 <- Clim_PCA$scores[,2]




#---------------------------------------------------------#
#### 4. Run Cox Ph model with fixed 10 day env windows ####
#---------------------------------------------------------#

## **NOTE** not used in the main paper ##


#--------------------------------------#
#### 4.1 Prepare data set for model ####
#--------------------------------------#

## remove rows with missing explanatories
Inc_br <- Inc_br %>% drop_na(staging_length)

## filter out the birds that deffered
Inc_br <- filter(Inc_br, is.na(attempt_start) == F)
Inc_br <- filter(Inc_br, Acc == "Y")

## examine which potential explanatories could be correlated
test2 <- subset(Inc_br, select = c("breeding_lat", "Green_centre", "laying_centre", "Comp1", "Comp2",
                                   "staging_length", "Gr10precip", "total_ODBA"))
correlation_matrix2 <- cor(test2)

## set explanatories as correct class
Inc_br$year <- as.factor(as.character(Inc_br$year))

## create nest fate column
Inc_br$Fate23 <- ifelse(Inc_br$success23 == 0, 1, 0)
Inc_br$Fate24 <- ifelse(Inc_br$success24 == 0, 1, 0)
Inc_br$high_lat <- ifelse(Inc_br$breeding_lat > 69.5, "high", "low")
Inc_br$high_lat <- as.factor(Inc_br$high_lat)


## scale the explanaotories
Inc_br_sc <- Inc_br %>% 
              mutate(Gr10precip = (Gr10precip-mean(Gr10precip))/sd(Gr10precip),
                     total_ODBA = (total_ODBA-mean(total_ODBA))/sd(total_ODBA),
                     staging_length = (staging_length-mean(staging_length))/sd(staging_length),
                     Green_centre = (Green_centre-mean(Green_centre))/sd(Green_centre),
                     mig_duration = (mig_duration-mean(mig_duration))/sd(mig_duration),
                     laying_centre = (laying_centre-mean(laying_centre))/sd(laying_centre),
                     failure_centre = (failure_centre-mean(failure_centre))/sd(failure_centre))


#-------------------------#
#### 4.2 Cox ph models ####
#-------------------------#

# time to event model
time_model24 <- coxph(Surv(length, Fate24) ~ Gr10precip + Gr10Sn + laying_centre + Green_centre + year,
                      data= Inc_br_sc)

summary(time_model24)
drop1(time_model24, test = "Chi")


## Checking Cox ph model assumptions 

## check deviance residuals vs predicted values
ggcoxdiagnostics(time_model24)
ggplot(data = NULL, aes(x=Inc_br_sc$year, y=resid(time_model24, type = "martingale"))) + geom_point() + geom_smooth()

## check for proportionality, using Schoenfeld residuals for each variable
ggcoxzph(cox.zph(time_model24)) #  For each covariate it produces plots with scaled Schoenfeld residuals against the time.


####---- Use MuMIn to perfrom model selection ----####

## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
ms2 <- MuMIn::dredge(time_model24, trace = 2)

## perform model selection
ms2_sub <- subset(ms2, !nested(.), recalc.weights=T)
ms2_sub <- subset(ms2_sub, delta <= 56, recalc.weights=T)
ms2_sub



#----------------------------------------#
#### 5. Counting process Cox ph model ####
#----------------------------------------#

#-----------------------------------#
#### 5.1 Prepare Incubation data ####
#-----------------------------------#

## create data set
Surv_tab <- as.data.frame(Inc_breeders$Tag_year)
Surv_tab$FirstFound <- yday(lubridate::dmy(Inc_breeders$attempt_start))
Surv_tab$LastChecked <- yday(lubridate::dmy(Inc_breeders$attempt_end))
Surv_tab$LastPresent <- as.numeric(yday(lubridate::dmy(Inc_breeders$attempt_end)))+1
Surv_tab$Fate <- ifelse(Inc_breeders$success23 == 0, 1, 0)
Surv_tab$AgeFound <- 1
Surv_tab$AgeDay1 <- Surv_tab$AgeFound - Surv_tab$FirstFound
Surv_tab$year <- year(lubridate::dmy(Inc_breeders$attempt_start))
Surv_tab$laying_centre <- Inc_breeders$laying_centre
Surv_tab$staging_length <- Inc_breeders$staging_length
Surv_tab$mig_duration <- Inc_breeders$mig_duration
Surv_tab$Green_centre <- Inc_breeders$Green_centre
Surv_tab$Acc <- Inc_breeders$Acc
Surv_tab$breeding_lat <- Inc_breeders$breeding_lat
Surv_tab$total_ODBA <- Inc_breeders$total_ODBA
Surv_tab$ID <- Inc_breeders$ID

## join on Env data windows from Greenland arrival
names(Surv_tab)[1] <- "Tag_year"
Surv_tag2 <- inner_join(Surv_tab, Window, by = "Tag_year")
Surv_tag2$Comp1 <- Clim_PCA$scores[,1] # add first two PCA axis
Surv_tag2$Comp2 <- Clim_PCA$scores[,2]
Surv_tag2 <- filter(Surv_tag2, is.na(FirstFound) == F)

## Use function from other script to expand out data set so each day of incubation has own row
Exp_ph <- expand.nest.data.ph(Surv_tag2)

## fix glitches in Fail and then recalculate the Surv column
Exp_ph$Fail <- ifelse(Exp_ph$Fate == 1 & (Exp_ph$LastPresent) == Exp_ph$End, 1, 0)
Exp_ph$Surv <- Surv(time=Exp_ph$Start, time2=Exp_ph$End, event=Exp_ph$Fail, type="counting")

## Now join on the rolling env windows for each day of the incubation
setnames(Exp_ph, old = "Start", new = "yday")
Exp_ph <- inner_join(Exp_ph, Roll, by = c("Tag_year", "yday"))

## remove rows with missing explanatories
Exp_ph <- Exp_ph %>% drop_na(c(staging_length))


## Create 2 day weather windows that move in time
## Prepare weather data; change which variable selected to change moving window size
Roll_sub <- subset(Roll, select=c("avg_temp5", "sum_precip5", "Tag_year", "yday"))
Exp_ph$cutoff <- ifelse(Exp_ph$breeding_lat> 69.5, "High", "low")
ggplot(Exp_ph, aes(x = yday, y = Tag_year, colour = sum_precip)) + geom_point() + facet_wrap(~year+cutoff)

# ## Creates staggered yday column
# Exp_ph$yday1 <- Exp_ph$yday-1; Exp_ph$yday2 <- Exp_ph$yday-2; Exp_ph$yday3 <- Exp_ph$yday-3; Exp_ph$yday4 <- Exp_ph$yday-4; Exp_ph$yday5 <- Exp_ph$yday-5; Exp_ph$yday6 <- Exp_ph$yday-6
# Exp_ph$yday7 <- Exp_ph$yday-7; Exp_ph$yday8 <- Exp_ph$yday-8; Exp_ph$yday9 <- Exp_ph$yday-9; Exp_ph$yday10 <- Exp_ph$yday-10; Exp_ph$yday11 <- Exp_ph$yday-11; Exp_ph$yday12 <- Exp_ph$yday-12
# Exp_ph$yday15 <- Exp_ph$yday-15
# ## Join weather data to staggered dates
# setnames(Roll_sub, old = c("yday", "avg_temp5", "sum_precip5"), new = c("yday1", "avg_temp2sub1", "sum_precip2sub1"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday1"))
# setnames(Roll_sub, old = c("yday1", "avg_temp2sub1", "sum_precip2sub1"), new = c("yday2", "avg_temp2sub2", "sum_precip2sub2"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday2"))
# setnames(Roll_sub, old = c("yday2", "avg_temp2sub2", "sum_precip2sub2"), new = c("yday3", "avg_temp2sub3", "sum_precip2sub3"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday3"))
# setnames(Roll_sub, old = c("yday3", "avg_temp2sub3", "sum_precip2sub3"), new = c("yday4", "avg_temp2sub4", "sum_precip2sub4"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday4"))
# setnames(Roll_sub, old = c("yday4", "avg_temp2sub4", "sum_precip2sub4"), new = c("yday5", "avg_temp2sub5", "sum_precip2sub5"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday5"))
# setnames(Roll_sub, old = c("yday5", "avg_temp2sub5", "sum_precip2sub5"), new = c("yday6", "avg_temp2sub6", "sum_precip2sub6"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday6"))
# setnames(Roll_sub, old = c("yday6", "avg_temp2sub6", "sum_precip2sub6"), new = c("yday7", "avg_temp2sub7", "sum_precip2sub7"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday7"))
# setnames(Roll_sub, old = c("yday7", "avg_temp2sub7", "sum_precip2sub7"), new = c("yday8", "avg_temp2sub8", "sum_precip2sub8"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday8"))
# setnames(Roll_sub, old = c("yday8", "avg_temp2sub8", "sum_precip2sub8"), new = c("yday9", "avg_temp2sub9", "sum_precip2sub9"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday9"))
# setnames(Roll_sub, old = c("yday9", "avg_temp2sub9", "sum_precip2sub9"), new = c("yday10", "avg_temp2sub10", "sum_precip2sub10"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday10"))
# setnames(Roll_sub, old = c("yday10", "avg_temp2sub10", "sum_precip2sub10"), new = c("yday11", "avg_temp2sub11", "sum_precip2sub11"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday11"))
# setnames(Roll_sub, old = c("yday11", "avg_temp2sub11", "sum_precip2sub11"), new = c("yday12", "avg_temp2sub12", "sum_precip2sub12"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday12"))
# setnames(Roll_sub, old = c("yday12", "avg_temp2sub12", "sum_precip2sub12"), new = c("yday15", "avg_temp2sub15", "sum_precip2sub15"))
# Exp_ph <- inner_join(Exp_ph, Roll_sub, by = c("Tag_year", "yday15"))

## set year as correct class
Exp_ph2 <- filter(Exp_ph, Acc == "Y")
Exp_ph2$year <- as.factor(Exp_ph2$year)
table(Exp_ph2$year)

## check correlation between various predictor variables
sub_cor <- subset(Exp_ph2, select = c("laying_centre", "Green_centre", "staging_length",
                                      "Comp1", "Comp2", "breeding_lat", "total_ODBA", "Gr10precip"))
co_matrix <- cor(sub_cor) ## NDVI windows starting from Greenland arrival and incubation end are correlated

## check for correlation between sub-population and continuous predictors
biserial.cor(Exp_ph2$Comp1, Exp_ph2$cutoff)
biserial.cor(Exp_ph2$Gr10precip, Exp_ph2$cutoff)
biserial.cor(Exp_ph2$laying_centre, Exp_ph2$cutoff)
biserial.cor(Exp_ph2$Green_centre, Exp_ph2$cutoff)
biserial.cor(Exp_ph2$avg_temp, Exp_ph2$cutoff)



#----------------------#
#### 5.2 Run models ####
#----------------------#

#### scale the explanatory
Exp_ph2_sc <- Exp_ph2 %>% 
              mutate(Gr10precip = (Gr10precip-mean(Gr10precip))/sd(Gr10precip),
                     total_ODBA = (total_ODBA-mean(total_ODBA))/sd(total_ODBA),
                     staging_length = (staging_length-mean(staging_length))/sd(staging_length),
                     Green_centre = (Green_centre-mean(Green_centre))/sd(Green_centre),
                     sum_precip = (sum_precip-mean(sum_precip))/sd(sum_precip),
                     avg_temp = (avg_temp-mean(avg_temp))/sd(avg_temp),
                     laying_centre = (laying_centre-mean(laying_centre))/sd(laying_centre),
                     Comp1 = (Comp1-mean(Comp1))/sd(Comp1),
                     sum_precip2 = (sum_precip2-mean(sum_precip2))/sd(sum_precip2),
                     avg_temp2 = (avg_temp2-mean(avg_temp2))/sd(avg_temp2),
                     sum_precip3 = (sum_precip3-mean(sum_precip3))/sd(sum_precip3),
                     avg_temp3 = (avg_temp3-mean(avg_temp3))/sd(avg_temp3),
                     sum_precip4 = (sum_precip4-mean(sum_precip4))/sd(sum_precip4),
                     avg_temp4 = (avg_temp4-mean(avg_temp4))/sd(avg_temp4),
                     sum_precip5 = (sum_precip5-mean(sum_precip5))/sd(sum_precip5),
                     avg_temp5 = (avg_temp5-mean(avg_temp5))/sd(avg_temp5),
                     sum_precip10 = (sum_precip10-mean(sum_precip10))/sd(sum_precip10),
                     avg_temp10 = (avg_temp10-mean(avg_temp10))/sd(avg_temp10))


#### Run models with varying  window length of time-dependent climatic variables 

mod.null <-  coxme(Surv ~ (1|ID), data= Exp_ph2)

mod_roll_ran  <-  coxme(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip + avg_temp + laying_centre + Green_centre + year + (1|ID), 
                        data= Exp_ph2_sc)
mod_roll_ran2  <-  coxme(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip2  + avg_temp2 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran3  <-  coxme(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip3  + avg_temp3 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran4  <-  coxme(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip4  + avg_temp4 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran5  <-  coxme(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip5  + avg_temp5 + laying_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran10  <-  coxme(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip10  + avg_temp10 + laying_centre + Green_centre + year + (1|ID), 
                          data= Exp_ph2_sc)

## compare the different models with AICc
AICc(mod.null, mod_roll_ran, mod_roll_ran2, mod_roll_ran3, mod_roll_ran4, mod_roll_ran5, mod_roll_ran10)

## get the summary from the best model
summary(mod_roll_ran)




#------------------------#
#### 5.3 Model Checks ####
#------------------------#

## check residuals vs predicted values
plot(predict(mod_roll_ran), resid(mod_roll_ran, type = "deviance"))

## can get the following model checks to work if treat individual ID as a frailty term instead of a random effect
mod_roll_check  <-  coxph(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip + avg_temp + laying_centre + Green_centre + year + frailty(ID), 
                        data= Exp_ph2_sc)
## check for proportionality, using Schoenfeld residuals for each variable
ggcoxzph(cox.zph(mod_roll_check)) #  For each covariate it produces plots with scaled Schoenfeld residuals against the time.



#---------------------------#
#### 5.3 Model Selection ####
#---------------------------#

## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
ms2 <- MuMIn::dredge(mod_roll_ran, trace = 2)

## apply nesting rule and select models within 6 AICc points of top model
ms2_sub <- subset(ms2, !nested(.), recalc.weights=T)
ms2_sub6 <- subset(ms2_sub, delta <= 6, recalc.weights=T)
ms2_sub6


## OBTAIN 95% CIs FOR TOP MODEL SET 

topmod1  <-  coxme(Surv ~ Gr10precip + cutoff + Green_centre + (1|ID), data= Exp_ph2_sc)
summary(topmod1);confint(topmod1)

topmod2  <-  coxme(Surv ~ cutoff + Green_centre + (1|ID), data= Exp_ph2_sc)
summary(topmod2);confint(topmod2)

topmod3 <-  coxme(Surv ~ cutoff + (1|ID), data= Exp_ph2_sc)
summary(topmod3);confint(topmod3)

topmod4 <-  coxme(Surv ~ Green_centre + Comp1 + avg_temp2 + (1|ID), data= Exp_ph2_sc)
summary(topmod4);confint(topmod4)

topmod5 <-  coxme(Surv ~ Green_centre + Comp1 + (1|ID), data= Exp_ph2_sc)
summary(topmod5);confint(topmod5)




#-----------------------------------------------------------#
#### 6. Create Kaplain-Meier for the two sub-populations ####
#-----------------------------------------------------------#

## First create the plot using ggadjustedcurves()
## Cant use mixed effects models or friality terms but survival rates are conditional on other terms in the model
## Doesn't seem to allow addition of 95% CIs

## Run the top model from the MUMin model selection
topmod_simp  <-  coxph(Surv ~ Gr10precip + cutoff + Green_centre, 
                          data= Exp_ph2_sc)

## use the ggadjustedcurves function to create the plot
survminer::ggadjustedcurves(topmod_simp, data = Exp_ph2_sc, variable = "cutoff", method = "conditional",
                            xlim = c(138, 192))



## Now create the plot using ggsurvplot()
## This allows you to add 95% CIs but isn't conditional on other variables in the model

## first create a survfit object with sub-popualtion as the only explanatory
Pop_fit <- survfit(Surv ~ cutoff, data = Exp_ph2_sc)

## Plot the two survival curves using the ggsurvplot function
ggsurvplot(Pop_fit, data = Exp_ph2_sc, conf.int = TRUE, pval = FALSE,
           legend = "right", legend.title = "Sub-population", legend.labs = c("Wexford", "Islay"),
           censor = F, conf.int.alpha = 0.15, xlab = "Day of year", 
           palette = c("#D55E00", "#0072B2"), break.x.by = 10,
           font.x = c(18, "black"), font.y = c(18, "black"), font.legend = c(14,"black"), xlim = c(144, 185))

## Save a plot
ggsave("Paper Plots/Figure 4- Kaplan Meir plot.png", 
       width = 22, height = 18, units = "cm")
