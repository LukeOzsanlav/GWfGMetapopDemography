## Luke Ozsanlav-Harris

## Updates: 19/08/2022

## Run Cox proportional hazard models on incubation data to model variaiton in nest failure rates
## Use climatic variables as the explanatories in the model and exmaine effect of sub-population
## Need to load the function from the script CoxPH-additional-functions
 

## packages required
pacman::p_load(tidyverse, data.table, survival, finalfit, ggpubr,
               survminer, coxme, MuMIn, forestplot, dplyr, ltm, rphylopic)

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
Roll$Tag_year <- Roll$tag_year
Roll$yday <- yday(as.Date(Roll$date))
Roll$ID <- NULL
Roll$year <- NULL
Roll$tag_year <- NULL




#---------------------------------------------#
#### 2. Read in Incubation data and format ####
#---------------------------------------------#

## Read in the incubation data set
Inc_ext <- fread("Outputs/Incubation_attempts_for_models.csv")
Inc_breeders <- Inc_ext

## add tag_date column with tag ID and incubation start date
Inc_breeders$tag_date <- paste0(Inc_breeders$ID, "_", Inc_breeders$attempt_end)

## filter out just the Ornitela data
Inc_breeders <- Inc_breeders %>% 
                mutate(tagtype = substr(ID, 1, 1)) %>% 
                filter(tagtype == 1) %>% 
                dplyr::select(-tagtype)


#-------------------------------------#
#### 3. PCA of climatic variables ####
#-------------------------------------#

## add 10 day env windows from greenland arrival to main data set
Inc_br <- inner_join(Inc_breeders, Window, by = "Tag_year") %>% drop_na(attempt_start)
Inc_br <- filter(Inc_br, !year == 2017)

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
# ggplot(loadings, aes(x=Symbol, y = Weight)) + 
#   geom_bar(stat='identity') + 
#   facet_grid(Component ~ ., scales = "free_y")

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

## filter out the birds that deferred
Inc_br <- filter(Inc_br, is.na(attempt_start) == F)

## examine which potential explanatories could be correlated
test2 <- subset(Inc_br, select = c("breeding_lat", "Green_centre", "laying_centre", "Comp1", "Comp2",
                                   "staging_length", "Gr10precip"))
correlation_matrix2 <- cor(test2)

## set explanatories as correct class and add extra columns for models
Inc_br <- Inc_br %>% 
          mutate(year = as.factor(as.character(year)),
                 Fate24 = ifelse(success24 == 0, 1, 0),
                 high_lat = ifelse(breeding_lat > 69.5, "high", "low"),
                 high_lat = as.factor(high_lat),
                 sub_pop = ifelse(Ringing.location == "WEXF" | Ringing.location == "SESK", "HVAN", "LOWLANDS"),
                 sub_pop = as.factor(sub_pop))

## Can just filter out the Islay and Wexford data at this stage if required
# Inc_br <- filter(Inc_br, Ringing.location %in% c("ISLA", "WEXF"))
# Inc_br$sub_pop <- as.factor(Inc_br$Ringing.location)


## scale the explanaotories
Inc_br_sc <- Inc_br %>% 
              mutate(Gr10precip = (Gr10precip-mean(Gr10precip, na.rm = T))/sd(Gr10precip, na.rm = T),
                     staging_length = (staging_length-mean(staging_length, na.rm = T))/sd(staging_length, na.rm = T),
                     Green_centre = (Green_centre-mean(Green_centre, na.rm = T))/sd(Green_centre, na.rm = T),
                     laying_centre = (laying_centre-mean(laying_centre, na.rm = T))/sd(laying_centre, na.rm = T),
                     failure_centre = (failure_centre-mean(failure_centre, na.rm = T))/sd(failure_centre, na.rm = T),
                     Comp1 = (Comp1-mean(Comp1, na.rm = T))/sd(Comp1, na.rm = T))



#----------------------------#
#### 4.2 Run Cox ph model ####
#----------------------------#

# time to event model
time_model24 <- coxme(Surv(length, Fate24) ~ Gr10precip + Comp1 + sub_pop + laying_centre + Green_centre + year + (1|ID),
                      data= Inc_br_sc)

summary(time_model24)
confint(time_model24)



#------------------------#
#### 4.3 Model checks ####
#------------------------#

## check deviance residuals vs predicted values
ggcoxdiagnostics(time_model24)
ggplot(data = NULL, aes(x=Inc_br_sc$year, y=resid(time_model24, type = "martingale"))) + geom_point() + geom_smooth()

## check for proportionality, using Schoenfeld residuals for each variable
ggcoxzph(cox.zph(time_model24)) #  For each covariate it produces plots with scaled Schoenfeld residuals against the time.




#---------------------------#
#### 4.4 Model selection ####
#---------------------------#

## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
ms2 <- MuMIn::dredge(time_model24, trace = 2)

## perform model selection
ms2_sub <- subset(ms2, !nested(.), recalc.weights=T)
ms2_sub <- subset(ms2_sub, delta <= 6, recalc.weights=T)
ms2_sub


## run the top models from the dredge to get the confidence intervals
time_model1 <- coxme(Surv(length, Fate24) ~ Green_centre + sub_pop + (1|ID),
                     data= Inc_br_sc)
confint(time_model1)
summary(time_model1)

time_model2 <- coxme(Surv(length, Fate24) ~ Green_centre + Comp1 + (1|ID),
                     data= Inc_br_sc)
confint(time_model2)
summary(time_model2)




#----------------------------------------#
#### 5. Counting process Cox ph model ####
#----------------------------------------#

#-----------------------------------#
#### 5.1 Prepare Incubation data ####
#-----------------------------------#

## create data set
Surv_tab <- as.data.frame(Inc_br$Tag_year)
Surv_tab$FirstFound <- yday(lubridate::ymd(Inc_br$attempt_start))
Surv_tab$LastChecked <- yday(lubridate::ymd(Inc_br$attempt_end))
Surv_tab$LastPresent <- as.numeric(yday(lubridate::ymd(Inc_br$attempt_end)))+1 
#*** USE this for nest age aligned analysis
# Surv_tab$FirstFound <- 0
# Surv_tab$LastChecked <- as.numeric(Inc_br$length)-1
# Surv_tab$LastPresent <- as.numeric(Inc_br$length)
#***
Surv_tab$Fate <- ifelse(Inc_br$success24 == 0, 1, 0)
Surv_tab$AgeFound <- 1
Surv_tab$AgeDay1 <- Surv_tab$AgeFound - Surv_tab$FirstFound
Surv_tab$year <- year(lubridate::ymd(Inc_br$attempt_start))
Surv_tab$laying_centre <- Inc_br$laying_centre
Surv_tab$staging_length <- Inc_br$staging_length
Surv_tab$mig_duration <- Inc_br$mig_duration
Surv_tab$Green_centre <- Inc_br$Green_centre
Surv_tab$Acc <- Inc_br$Acc
Surv_tab$breeding_lat <- Inc_br$breeding_lat
Surv_tab$total_ODBA <- Inc_br$total_ODBA
Surv_tab$ID <- Inc_br$ID
Surv_tab$RingLoc <- Inc_br$Ringing.location
Surv_tab$attempt_start <- Inc_br$attempt_start

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
#Exp_ph$yday <- yday(lubridate::ymd(Exp_ph$attempt_start) + lubridate::days(Exp_ph$NestAge)) #*** USE this for nest age aligned analysis
Exp_ph <- Exp_ph %>%  dplyr::rename(yday = Start)
Exp_ph2 <- inner_join(Exp_ph, Roll, by = c("Tag_year", "yday"))

## Create cut off variable
Exp_ph2$cutoff <- ifelse(Exp_ph2$breeding_lat> 69.5, "High", "low")

## set year as correct class
Exp_ph2$year <- as.factor(Exp_ph2$year)
table(Exp_ph2$year)

## check correlation between various predictor variables
sub_cor <- subset(Exp_ph2, select = c("laying_centre", "Green_centre", "staging_length",
                                      "Comp1", "Comp2", "breeding_lat", "Gr10precip"))
co_matrix <- cor(sub_cor) ## NDVI windows starting from Greenland arrival and incubation end are correlated

## Add on a sub-population variable
Exp_ph2$sub_pop <- ifelse(Exp_ph2$RingLoc == "WEXF" | Exp_ph2$RingLoc == "SESK", "HVAN", "LOWLANDS")
Exp_ph2$sub_pop <- as.factor(Exp_ph2$sub_pop)

## check for correlation between sub-population and continuous predictors
biserial.cor(Exp_ph2$Comp1, Exp_ph2$sub_pop)
biserial.cor(Exp_ph2$Gr10precip, Exp_ph2$sub_pop)
biserial.cor(Exp_ph2$laying_centre, Exp_ph2$sub_pop)
biserial.cor(Exp_ph2$Green_centre, Exp_ph2$sub_pop)
biserial.cor(Exp_ph2$avg_temp, Exp_ph2$sub_pop)
biserial.cor(Exp_ph2$sum_precip, Exp_ph2$sub_pop)

ggplot() + geom_boxplot(data = Exp_ph2, aes(x= sub_pop, y =sum_precip)) + theme_bw()
hist(Exp_ph2$sum_precip)



#----------------------#
#### 5.2 Run models ####
#----------------------#

#### scale the explanatory
Exp_ph2_sc <- Exp_ph2 %>% 
              mutate(Gr10precip = (Gr10precip-mean(Gr10precip))/sd(Gr10precip),
                     staging_length = (staging_length-mean(staging_length, na.rm=T))/sd(staging_length, na.rm=T),
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
Exp_ph2_sc$year <- as.character(Exp_ph2_sc$year)

mod.null <-  coxme(Surv ~ (1|ID), data= Exp_ph2_sc)

mod_roll_ran  <-  coxme(Surv ~ Gr10precip + Comp1 + sub_pop + sum_precip + avg_temp + laying_centre + Green_centre + year + (1|ID), 
                        data= Exp_ph2_sc)
mod_roll_ran2  <-  coxme(Surv ~ Gr10precip + Comp1 + sub_pop + sum_precip2 + avg_temp2 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran3  <-  coxme(Surv ~ Gr10precip + Comp1 + sub_pop + sum_precip3 + avg_temp3 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran4  <-  coxme(Surv ~ Gr10precip + Comp1 + sub_pop + sum_precip4 + avg_temp4 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran5  <-  coxme(Surv ~ Gr10precip + Comp1 + sub_pop + sum_precip5 + avg_temp5 + laying_centre + + Green_centre + year + (1|ID), 
                         data= Exp_ph2_sc)
mod_roll_ran10  <-  coxme(Surv ~ Gr10precip + Comp1 + sub_pop + sum_precip10 + avg_temp10 + laying_centre + Green_centre + year + (1|ID), 
                          data= Exp_ph2_sc)

## compare the different models with AICc
AICc(mod.null, mod_roll_ran, mod_roll_ran2, mod_roll_ran3, mod_roll_ran4, mod_roll_ran5, mod_roll_ran10)

## get the summary from the best model
summary(mod_roll_ran)
confint(mod_roll_ran)

## try and visualize why we might be getting this precip effect
# summary(Exp_ph2_sc$sum_precip)
# ggplot(Exp_ph2_sc) + geom_point(aes(y=Tag_year, x= yday, colour = sum_precip)) + theme_bw() + scale_color_viridis_b()
# hist(Exp_ph2$sum_precip)
# Precip_check <- filter(Exp_ph2, sum_precip > 0.1)
# biserial.cor(Precip_check$sum_precip, Precip_check$sub_pop)


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
#### 5.4 Model Selection ####
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

topmod1  <-  coxme(Surv ~ avg_temp + Green_centre + sum_precip + (1|ID), data= Exp_ph2_sc)
summary(topmod1);confint(topmod1)

topmod2  <-  coxme(Surv ~ Green_centre + sub_pop + sum_precip + (1|ID), data= Exp_ph2_sc)
summary(topmod2);confint(topmod2)

topmod3 <-  coxme(Surv ~ Comp1 + Green_centre + sum_precip + (1|ID), data= Exp_ph2_sc)
summary(topmod3);confint(topmod3)

topmod4 <-  coxme(Surv ~ sub_pop + sum_precip + + (1|ID), data= Exp_ph2_sc)
summary(topmod4);confint(topmod4)

topmod5 <-  coxme(Surv ~ Green_centre + sum_precip + (1|ID), data= Exp_ph2_sc)
summary(topmod5);confint(topmod5)

topmod6  <-  coxme(Surv ~ avg_temp + sum_precip + (1|ID), data= Exp_ph2_sc)
summary(topmod6);confint(topmod6)




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
Pop_fit <- survfit(Surv ~ sub_pop, data = Exp_ph2_sc)

## Plot the two survival curves using the ggsurvplot function
ggsurvplot(Pop_fit, data = Exp_ph2_sc, conf.int = TRUE, pval = FALSE,
           legend = "right", legend.title = "Sub-population", legend.labs = c("Wexford", "Islay"),
           censor = F, conf.int.alpha = 0.15, xlab = "Day of year", 
           palette = c("#D55E00", "#0072B2"), break.x.by = 10,
           font.x = c(18, "black"), font.y = c(18, "black"), font.legend = c(14,"black"), xlim = c(144, 185))

## Save a plot
ggsave("Paper Plots/Figure 5- Kaplan Meir plot.png", 
       width = 22, height = 18, units = "cm")









#-------------------------------------#
#### 7. Run populations separately ####
#-------------------------------------#

#------------------------------------------#
#### 7.1 Filter out the two populations ####
#------------------------------------------#

## filter the data into the two sub-populations
Exp_ph2_Hvan <- filter(Exp_ph2_sc, sub_pop == "HVAN")
Exp_ph2_Low <- filter(Exp_ph2_sc, sub_pop == "LOWLANDS")

# how manny bird years in each analysis
length(unique(Exp_ph2_Hvan$Tag_year))
length(unique(Exp_ph2_Low$Tag_year))

## average temp during incubation




#---------------------------#
#### 7.2 Run Wexf models ####
#---------------------------#


mod.nullWexf <-  coxme(Surv ~ (1|ID), data= Exp_ph2_Hvan)

mod_Wexf_ran  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip + avg_temp + laying_centre + Green_centre + year + (1|ID), 
                        data= Exp_ph2_Hvan)
mod_Wexf_ran2  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip2 + avg_temp2 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Hvan)
mod_Wexf_ran3  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip3 + avg_temp3 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Hvan)
mod_Wexf_ran4  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip4 + avg_temp4 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Hvan)
mod_Wexf_ran5  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip5 + avg_temp5 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Hvan)
mod_Wexf_ran10  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip10 + avg_temp10 + laying_centre + Green_centre + year + (1|ID), 
                          data= Exp_ph2_Hvan)

## compare the different models with AICc
AICc(mod.nullWexf, mod_Wexf_ran, mod_Wexf_ran2, mod_Wexf_ran3, mod_Wexf_ran4, mod_Wexf_ran5, mod_Wexf_ran10)

## get the summary from the best model
summary(mod_Wexf_ran4)
confint(mod_Wexf_ran4)


#-----------------------------#
#### 7.3 Wexf Model Checks ####
#-----------------------------#

## check residuals vs predicted values
plot(predict(mod_Wexf_ran), resid(mod_Wexf_ran, type = "deviance"))

## can get the following model checks to work if treat individual ID as a frailty term instead of a random effect
mod_Wexf_check  <-  coxph(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip + avg_temp + laying_centre + Green_centre + year + frailty(ID), 
                          data= Exp_ph2_Hvan)
## check for proportionality, using Schoenfeld residuals for each variable
ggcoxzph(cox.zph(mod_Wexf_check)) #  For each covariate it produces plots with scaled Schoenfeld residuals against the time.



#-------------------------------#
#### 7.4 Wexf Model Plotting ####
#-------------------------------#


## get the estimates and confidence intervals from the models
Ests <- data.frame(summary(mod_Wexf_ran4)$coefficients)
Ests <- Ests %>% cbind(confint(mod_Wexf_ran4)) %>% dplyr::rename(Estimate = `summary.mod_Wexf_ran4..coefficients`)
rownames(Ests) <- c("ArrPrecip", "Clim", "Precip(t)", "Temp(t)",
                    "Inc", "Arrival", "year[2019]", "year[2020]", "year[2021]")

## make a forest plot fo the model estimates
Ests <- Ests %>% 
        map_df(rev) %>% 
        mutate(Param = rev(rownames(Ests)),
               Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))

## get goose image from rphylopic
goose <- name_search(text = "Anser albifrons", options = "namebankID")[[1]]
goose_img <- image_data(name_images(uuid = goose$uid[1])$supertaxa[[1]]$uid, size=1024)[[1]]

Wexf <- ggplot(Ests) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
        geom_errorbarh(aes(y=Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, size =0.5) +
        geom_point(aes(y=Param, x= Estimate, color = Sig), size = 2.5) +
        theme_bw() +
        ylab("") +
        xlim(-5, 5) +
        scale_color_manual(values=c("#B2BABB", "#E74C3C")) +
        annotate(geom="text", x=4.3, y=9.2, label="Wexford",color="#D55E00", size =7) +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.title=element_text(size=16), 
              legend.title=element_text(size=14),
              axis.text=element_text(size=14), 
              legend.text=element_text(size=12),
              panel.grid.minor.x = element_blank(),
              legend.position = "none") +
        add_phylopic(goose_img, x=-4.5, y=2, ysize = 1.6, alpha=1, color="#D55E00")





#----------------------------#
#### 7.5 Run Islay Models ####
#----------------------------#

mod.null <-  coxme(Surv ~ (1|ID), data= Exp_ph2_Low)

mod_Islay_ran  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip + avg_temp + laying_centre + Green_centre + year + (1|ID), 
                        data= Exp_ph2_Low)
mod_Islay_ran2  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip2 + avg_temp2 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Low)
mod_Islay_ran3  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip3 + avg_temp3 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Low)
mod_Islay_ran4  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip4 + avg_temp4 + laying_centre + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Low)
mod_Islay_ran5  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip5 + avg_temp5 + laying_centre + + Green_centre + year + (1|ID), 
                         data= Exp_ph2_Low)
mod_Islay_ran10  <-  coxme(Surv ~ Gr10precip + Comp1 + sum_precip10 + avg_temp10 + laying_centre + Green_centre + year + (1|ID), 
                          data= Exp_ph2_Low)

## compare the different models with AICc
AICc(mod.null, mod_Islay_ran, mod_Islay_ran2, mod_Islay_ran3, mod_Islay_ran4, mod_Islay_ran5, mod_Islay_ran10)

## get the summary from the best model
summary(mod_Islay_ran)
confint(mod_Islay_ran)


#------------------------------#
#### 7.6 Islay Model Checks ####
#------------------------------#

## check residuals vs predicted values
plot(predict(mod_Islay_ran), resid(mod_Islay_ran, type = "deviance"))

## can get the following model checks to work if treat individual ID as a frailty term instead of a random effect
mod_Islay_check  <-  coxph(Surv ~ Gr10precip + Comp1 + cutoff + sum_precip + avg_temp + laying_centre + Green_centre + year + frailty(ID), 
                          data= Exp_ph2_Low)
## check for proportionality, using Schoenfeld residuals for each variable
ggcoxzph(cox.zph(mod_Islay_check)) #  For each covariate it produces plots with scaled Schoenfeld residuals against the time.



#--------------------------------#
#### 7.7 Islay Model Plotting ####
#--------------------------------#


## get the estimates and confidence intervals from the models
Ests2 <- data.frame(summary(mod_Islay_ran)$coefficients)
Ests2 <- Ests2 %>% cbind(confint(mod_Islay_ran)) %>% dplyr::rename(Estimate = `summary.mod_Islay_ran..coefficients`)
rownames(Ests2) <- c("ArrPrecip", "Clim", "Precip(t)", "Temp(t)",
                    "Inc", "Arrival", "year[2020]", "year[2021]")

## make a forest plot fo the model estimates
Ests2 <- Ests2 %>% 
  map_df(rev) %>% 
  mutate(Param = rev(rownames(Ests2)),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))

Islay <- ggplot(Ests2) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
          geom_errorbarh(aes(y=Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, size =0.5) +
          geom_point(aes(y=Param, x= Estimate, color = Sig), size = 2.5) +
          theme_bw() +
          ylab("") +
          xlim(-2.5, 2.5) +
          annotate(geom="text", x=2.25, y=8.2, label="Islay",color="#0072B2", size =7) +
          scale_color_manual(values=c("#B2BABB", "#E74C3C")) +
          theme(panel.grid.minor.y = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_blank(),
                axis.title=element_text(size=16), 
                legend.title=element_text(size=14),
                axis.text=element_text(size=14), 
                legend.text=element_text(size=12),
                panel.grid.minor.x = element_blank(),
                legend.position = "none") +
          add_phylopic(goose_img, x=-2.2, y=2, ysize = 1.5, alpha=1, color="#0072B2")


## save the plot as a png
ggarrange(Wexf, Islay, ncol=1, nrow=2)
ggsave("Paper Plots/Figure 4- Nest survival forest plot.png",
       width = 22, height = 28, units = "cm")






#-------------------------------------#
#### 7.8  Model Sub-pop difference ####
#-------------------------------------#


## Run model with ust sub-populaiton in
mod_roll_ran  <-  coxme(Surv ~ sub_pop + year + (1|ID), 
                        data= Exp_ph2_sc)

## create summary of model and CIs
summary(mod_roll_ran)
confint(mod_roll_ran)


