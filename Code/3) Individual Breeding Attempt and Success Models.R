## Luke Ozsanlav-Harris

## Run binomial GLMs to determine if any factors predict whether female GWfG
## attempt/defer breeding and if they do attempt, whether they are successful or fail


## Created: 26/08/2020

####        Script Contents       ####
## 1. Read in Averaged Env data
## 2. Read in Incubation data and format
## 3. Prepare data set for attempt/defer model
## 4. Binomal GLM for attempt/defer of breeding
## 5. Prepare dat set for success fail model
## 6. Binomal GLM for success/fail of breeding


## Changes that made: 01/10/2020
##  1. Use MuMIn to obtain more accurate p-values and conduct model simplification
##  2. Use a PCA to combine multiple climatic variables **DONE**
##  3. Include breeding latitude in models ** now breeding lat correlated with COMP1 and COMP2 **



## packages required
library(tidyverse)
library(data.table)
library(lme4)
library(DHARMa)
library(MuMIn)





######## 
## 1. ##
######## Read in Averaged Env data ####

setwd("~/PhD Documents/1_Tracking Data Chapters/Incuabtion Lengths/Combining Env data/Breeding Env data windows")
Window <- fread("Average_Env_data_from_arrival_windows_equalfixes.csv")

setwd("~/PhD Documents/1_Tracking Data Chapters/Incuabtion Lengths/Combining Env data/Daily rolling Env data")
Roll <- fread("Rolling_Env_data_per_day_breeding_season_evenfixes.csv")

setwd("~/PhD Documents/1_Tracking Data Chapters/Incuabtion Lengths/Combining Env data/Day centred env data")
Window_fail <- fread("Failure_env_windows.csv")

## Add extra columns to Roll data set so will bind later on to incubation data
Roll <- tidyr::separate(Roll, tag_date, into = c("tag", "date"), sep = "_")
Roll$Tag_year <- paste0(Roll$ID, "_", year(as.Date(Roll$date)))
Roll$yday <- yday(as.Date(Roll$date))





########
## 2. ##
######## Read in Incubation data and format ####

## Read in the incubation data set
setwd("~/PhD Documents/1_Tracking Data Chapters/Incuabtion Lengths/Breeding attempts")
Inc_ext <- fread("Incubation_attempts_extra.csv")
Inc_ext <- filter(Inc_ext, !Tag_year == "WHIT01_2018")

## filter out Acc tags only
Inc_acc <- filter(Inc_ext, Acc == "Y")

## add the Env data windows to the incubation data
Window$Tag_year <- as.character(Window$Tag_year)
Inc_env <- inner_join(Inc_ext, Window, by = "Tag_year")
Inc_env <- full_join(Inc_env, Window_fail, by = "Tag_year")






########
## 3. ##
######## PCA of climatic variables ####

## examine correlation between climatic variables first
sub1 <- subset(Inc_env, select = c("Gr10precip", "Gr10NDVI", "Gr10Sn", "Gr10temp", "Gr20precip", "Gr20NDVI", "Gr20Sn", "Gr20temp"))
cor_matrix <- cor(sub1)
## Preip doesnt correlate with any other climatic variable, absolute value of r < 0.25
## NDVI and Snow cover have high correlation, absolute  r > 0.75
## Temp has mid correlation with NDVI and snow cover, absolute r c0.5


## Select the varibale for PCA and scle them
sub_pca <- subset(Inc_env, select = c("Gr10NDVI", "Gr10Sn", "Gr10temp"))
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

## First two components explain 95% of the varation, will try and add both of them to binary models

## Add first and second principal component to data set
Inc_env$Comp1 <- Clim_PCA$scores[,1]
Inc_env$Comp2 <- Clim_PCA$scores[,2]


## Add 10 day precip from end
Roll_sub <- subset(Roll, select = c("Tag_year", "sum_precip10", "yday"))
Inc_env$yday <- yday(as.Date(Inc_env$attempt_end))
Inc_env2 <- inner_join(Inc_env, Roll_sub, by = c("Tag_year", "yday"))






########
## 4. ##
######## Binomal GLM for attempt/defer of breeding

## re-scale the explanatories
defer_sc <- Inc_env # rename man data set for later use
defer_sc[, c("Gr10precip", "total_ODBA", "staging_length", "Green_centre")] <- scale(
  defer_sc[, c("Gr10precip", "total_ODBA", "staging_length", "Green_centre")])

## remove rows with missing explanatories
defer_sc <- defer_sc %>% drop_na(c(total_ODBA, staging_length)) # remove two indivdiauls due to lack of migration duration

## check for correlations between possible explantory variables
test1 <- subset(defer_sc, select = c("Gr10precip", "Comp1", "Comp2", "Green_centre", "total_ODBA", "staging_length"))
correlation_matrix1 <- cor(test1) # can't have Green centre and Green yday, Green yday and Green10NDVI have cor of 0.47

## Make sure explanatories correct object class
defer_sc$year <- as.factor(defer_sc$year)

defer_sc$sub_pop <- ifelse(defer_sc$Ringing.location == "WEXF" | defer_sc$Ringing.location == "SESK", "HVAN", "LOWLANDS")


## run the binomial GLM
def_model <- glm(attempt ~  Gr10precip + Comp1 + Green_centre + total_ODBA + staging_length + year,
                 data = defer_sc,
                 family = binomial(link = "logit"))

## Can't get the random effects model to run with staging length in it
def_model_ran <- glmer(attempt ~  Gr10precip + Comp1 + Green_centre + sub_pop + year + (1|ID),
                 data = defer_sc,
                 family = binomial(link = "logit"),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
AICc(def_model)
AICc(def_model_ran)
## examine model
summary(def_model_ran)
drop1(def_model_ran, test = "Chi") # Liklihood ratio test, no terms significant


## test model assumptions using DHARMa
## simulate randomized quantile residuals
simulationOutput <- simulateResiduals(fittedModel = def_model, plot = T) # not violated

## plot the residuals against other predictors
par(mfrow=c(3,2))
plotResiduals(simulationOutput, form = defer_sc$Gr10precip) # not violated
plotResiduals(simulationOutput, form = defer_sc$mig_duration) # not violated
plotResiduals(simulationOutput, form = defer_sc$staging_length) # not violated
plotResiduals(simulationOutput, form = defer_sc$Green_centre) # not violated
plotResiduals(simulationOutput, form = defer_sc$Comp1) # not violated
plotResiduals(simulationOutput, form = defer_sc$omp2) # not violated

## test for outliers, dispersion and confomrity to expected distribution
testResiduals(simulationOutput) # not violated


#### Use MuMIn to perfrom model selection ####

## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
msdef <- MuMIn::dredge(def_model_ran, trace = 2)

# Visualize the model selection table: provides an idea of which variables in top ranked models
par(mfrow=c(1,1)); par(mar = c(3,5,6,4))
plot(msdef, labAsExpr = TRUE)
importance(msdef) # value of important of the different predictors

## perform model averaging
# select models within 5 AICc points of top model
msdef_sub <- subset(msdef, !nested(.), recalc.weights=T)
msdef_sub <- subset(msdef_sub, delta <= 5, recalc.weights=T)
avg_mods2 <- model.avg(msdef_sub, revised.var=T)

## Summarise outputs of model averaging
summary(avg_mods2)
confint(avg_mods2)






########
## 5. ##
######## Binomal GLM for succes/fail of attempted breeding

## filter out birds that did not attmept to breed
breeders <- dplyr::filter(Inc_env2, is.na(attempt_start) == F)

## set response as correct class 
breeders$success23 <- as.factor(breeders$success23)
breeders$success24 <- as.factor(breeders$success24)
breeders$year <- as.factor(breeders$year)

## scale the explanatories
breed_sc <- breeders
breed_sc[, c("Gr10precip", "breeding_lat", "total_ODBA", "staging_length", "Green_centre", "laying_centre", "sum_precip10", "fail5temp", "fail10temp", "fail15temp",       
             "fail5precip","fail10precip","fail15precip")] <- scale(
  breed_sc[, c("Gr10precip", "breeding_lat", "total_ODBA", "staging_length", "Green_centre", "laying_centre", "sum_precip10", "fail5temp", "fail10temp", "fail15temp",       
               "fail5precip","fail10precip","fail15precip")])

## remove rows with missing explanatories
breed_sc <- breed_sc %>% drop_na(c(staging_length, total_ODBA)) # lose two individuals

## check for correlations between possible explantory variables
test2 <- subset(breed_sc, select = c("Gr10precip", "breeding_lat", "total_ODBA", "staging_length", 
                                     "Green_centre", "laying_centre", "sum_precip10", "Comp1", "Comp2", "fail5temp", "fail10temp", "fail15temp",       
                                     "fail5precip", "fail10precip", "fail15precip"))
correlation_matrix2 <- cor(test2) # Note: breedng lat has moderate negative correlation with two PCA axis





####---- 23 day model ----####

# ## run the binimal GLM with 23 days as success cut off
# breed_mod23 <- glm(success23 ~ Gr10precip + Comp1 + laying_centre + Green_centre + 
#                                total_ODBA + staging_length + sum_precip10 + year,
#                              data = breed_sc,
#                              family = binomial(link = "logit"))
# summary(breed_mod23)
# drop1(breed_mod23, test = "Chi") #liklihood ratio test
#  
# ## check model assumptions with DHARMa
# ## simulate randomized quantile residuals
# simulationOutput23 <- simulateResiduals(fittedModel = breed_mod23, plot = T) # QQ plot not violated, lower quantile in second plot looks weird. Adding staging length^2 doesnt make plot look better buit stops signficant violation below
# 
# ## plot the residuals against other predictors
# par(mfrow=c(4,2))
# plotResiduals(simulationOutput23, form = breed_sc$Gr10precip) # not violated
# plotResiduals(simulationOutput23, form = breed_sc$Green_centre) # not violated shows slight quadratic shape in plot but probaly fine to leave
# plotResiduals(simulationOutput23, form = breed_sc$laying_centre) # not violated
# plotResiduals(simulationOutput23, form = breed_sc$mig_duration) # not violated
# plotResiduals(simulationOutput23, form = breed_sc$staging_length) # violation significant in orignal model, added in staging length^2 and sorted probelm
# plotResiduals(simulationOutput23, form = breed_sc$Comp1)
# plotResiduals(simulationOutput23, form = breed_sc$Comp2)
# 
# ## test for outliers, dispersion and confomrity to expected distribution
# testResiduals(simulationOutput23) # not violated
# 
# 
# #### Use MuMIn to perfrom model selection ####
# 
# ## change default "na.omit" to prevent models being fitted to different datasets
# options(na.action = "na.fail") 
# 
# ## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
# ms1 <- MuMIn::dredge(breed_mod23, trace = 2)
# 
# # Visualize the model selection table: provides an idea of which varibales in top ranked models
# par(mfrow=c(1,1)); par(mar = c(3,5,6,4))
# plot(ms1, labAsExpr = TRUE)
# importance(ms1) # value of important of the different predictors
# 
# ## perform model averaging
# # select models within 5 AICc points of top model
# ms1_sub <- subset(ms1, delta <= 5, recalc.weights=T)
# avg_mods <- model.avg(ms1_sub, revised.var=T)
# 
# ## Summarise outputs of model averaging
# summary(avg_mods)
# confint(avg_mods)





####---- 24 day model ----####
colnames(breed_sc)
table(breed_sc$year)
breed_sc$sub_pop <- ifelse(breed_sc$Ringing.location == "WEXF" | breed_sc$Ringing.location == "SESK", "HVAN", "LOWLANDS")

## run the binomial GLM with 24 days as success cut off
breed_mod24 <- glm(success24 ~ Comp1 + Gr10precip + laying_centre + Green_centre + staging_length ,
                   data = breed_sc,
                   family = binomial(link = "logit"))
breed_mod24_ran <- glmer(success24 ~ Comp1 + Gr10precip + laying_centre + Green_centre + sub_pop + year + (1|ID),
                   data = breed_sc,
                   family = binomial(link = "logit"),
                   control=glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=2e5)))
AICc(breed_mod24); AICc(breed_mod24_ran)

summary(breed_mod24_ran)
drop1(breed_mod24, test = "Chi") #liklihood ratio test


## check model assumptions with DHARMa
## simulate randomized quantile residuals
simulationOutput24 <- simulateResiduals(fittedModel = breed_mod24_ran, plot = T) # QQ plot not violated, lower quantile in second plot looks weird. Adding staging length^2 made plot look better

## plot the residuals against other predictors
par(mfrow=c(4,2))
plotResiduals(simulationOutput24, form = breed_sc$Gr10precip) # not violated
plotResiduals(simulationOutput24, form = breed_sc$breeding_lat) # not violated
plotResiduals(simulationOutput24, form = breed_sc$Green_centre) # not violated shows slight quadratic shape in plot but probaly fine to leave
plotResiduals(simulationOutput24, form = breed_sc$laying_centre) # not violated
plotResiduals(simulationOutput24, form = breed_sc$mig_duration) # not violated
plotResiduals(simulationOutput24, form = breed_sc$staging_length) # violation significant in orignal model, added in staging length^2 and sorted probelm
plotResiduals(simulationOutput24, form = breed_sc$Comp1)
plotResiduals(simulationOutput24, form = breed_sc$Comp2)

## test for outliers, dispersion and confomrity to expected distribution
testResiduals(simulationOutput24) # not violated


#### Use MuMIn to perfrom model selection ####

## change default "na.omit" to prevent models being fitted to different datasets
options(na.action = "na.fail") 

## create all candidate models using dredge, specify any dependencies (trace shows progress bar)
ms2 <- MuMIn::dredge(breed_mod24_ran, trace = 2)

# Visualize the model selection table: provides an idea of which variables in top ranked models
par(mfrow=c(1,1)); par(mar = c(3,5,6,4))
plot(ms2, labAsExpr = TRUE)
importance(ms2) # value of important of the different predictors

## perform model averaging
# select models within 5 AICc points of top model
ms2_sub <- subset(ms2, !nested(.), recalc.weights=T)
ms2_sub <- subset(ms2_sub, delta <= 6, recalc.weights=T)
avg_mods2 <- model.avg(ms2_sub, revised.var=T)

## Summarise outputs of model averaging
summary(avg_mods2)
confint(avg_mods2)







## Extract just the top model from AICc dredge
top_model24 <- get.models(ms2, subset = 1)[[1]]
summary(top_model24)



## run models that remained in the top models sets after model selection
## The extract the 95% CIs for the parameter estimates
mod1 <- glmer(success24 ~ sub_pop + (1|ID),
              data = breed_sc,
              family = binomial(link = "logit"),
              control=glmerControl(optimizer="bobyqa",
                                   optCtrl=list(maxfun=2e5)))
confint(mod1)
summary(mod1)


mod2 <- glmer(success24 ~ Green_centre + (1|ID),
              data = breed_sc,
              family = binomial(link = "logit"),
              control=glmerControl(optimizer="bobyqa",
                                   optCtrl=list(maxfun=2e5)))
confint(mod2)
