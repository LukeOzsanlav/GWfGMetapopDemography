## Luke Ozsanlav-Harris

## Run binomial GLMs to determine if any climatic factors and phenolgoy predict breeding propensity and success in female GWfG

## Created: 26/08/2020
## Updated: 19/08/2022
## Updated for review: 23/03/2023

## loadpackages required
pacman::p_load(tidyverse, data.table, lme4, DHARMa, MuMIn, glmmTMB, ltm)


#------------------------------------#
#### 1. Read in Averaged Env data ####
#------------------------------------#

## read in data sets
Window <- fread("Outputs/Average_Env_data_from_arrival_windows_equalfixes.csv")
Roll <- fread("Outputs/Rolling_Env_data_per_day_breeding_season_evenfixes.csv")

## Add extra columns to Roll data set so will bind later on to incubation data
Roll <- tidyr::separate(Roll, tag_date, into = c("tag", "date"), sep = "_")
Roll$Tag_year <- paste0(Roll$ID, "_", year(as.Date(Roll$date)))
Roll$yday <- yday(as.Date(Roll$date))




#---------------------------------------------#
#### 2. Read in Incubation data and format ####
#---------------------------------------------#

## Read in the incubation data set
Inc_ext <- fread("Outputs/Incubation_attempts_for_models.csv")

## add the Env data windows to the incubation data
Window$Tag_year <- as.character(Window$Tag_year)
Inc_acc <- inner_join(Inc_ext, Window, by = "Tag_year")

## filter for the years with Env data, these are the ones used in the original paper draft
## If i updated the env data then the data set would be larger
Inc_acc <- Inc_acc %>%  drop_na(Gr10NDVI, Gr10Sn, Gr10temp)

## drop the Ecotone tags
Inc_acc <- Inc_acc %>% 
           mutate(tagtype = substr(ID, 1, 1)) %>% 
           filter(tagtype == 1)





#------------------------------------#
#### 3. PCA of climatic variables ####
#------------------------------------#

## examine correlation between climatic variables first
sub1 <- subset(Inc_acc, select = c("Gr10precip", "Gr10NDVI", "Gr10Sn", "Gr10temp", "Gr20precip", "Gr20NDVI", "Gr20Sn", "Gr20temp"))
cor_matrix <- cor(sub1)
## Preip doesnt correlate with any other climatic variable, absolute value of r < 0.25
## NDVI and Snow cover have high correlation, absolute  r > 0.75
## Temp has mid correlation with NDVI and snow cover, absolute r c0.5


## Select the varibale for PCA and scle them
sub_pca <- subset(Inc_acc, select = c("Gr10NDVI", "Gr10Sn", "Gr10temp"))
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


## Add first and second principal component to data set
Inc_acc$Comp1 <- Clim_PCA$scores[,1]
Inc_acc$Comp2 <- Clim_PCA$scores[,2]


## Add 10 day precip from end
Roll_sub <- subset(Roll, select = c("Tag_year", "sum_precip10", "yday"))
Inc_acc$yday <- yday(as.Date(Inc_acc$attempt_end))
Inc_acc2 <- inner_join(Inc_acc, Roll_sub, by = c("Tag_year", "yday"))




#-----------------------------------------------#
#### 4. Binomal GLM for attempt/defer of breeding
#-----------------------------------------------#

#---------------------------#
#### 4.1 Format the data ####
#---------------------------#

## re-scale the explanatories
defer_sc <- Inc_acc # rename main data set for later use
defer_sc <- defer_sc %>% 
            mutate(Gr10precip = (Gr10precip-mean(Gr10precip, na.rm=T))/sd(Gr10precip, na.rm=T),
                   staging_length = (staging_length-mean(staging_length, na.rm=T))/sd(staging_length, na.rm=T),
                   Green_centre = (Green_centre-mean(Green_centre, na.rm=T))/sd(Green_centre, na.rm=T))


## check for correlations between possible explantory variables
test1 <- subset(defer_sc, select = c("Gr10precip", "Comp1", "Comp2", "Green_centre", "staging_length"))
correlation_matrix1 <- cor(test1) # can't have Green centre and Green yday, Green yday and Green10NDVI have cor of 0.47

## Make sure explanatories correct object class
defer_sc$year <- as.factor(defer_sc$year)

## define population column
defer_sc$sub_pop <- ifelse(defer_sc$Ringing.location == "WEXF" | defer_sc$Ringing.location == "SESK", "HVAN", "LOWLANDS")


#------------------------#
#### 4.2 Data Summary ####
#------------------------#

## get the summary of the data used in these models
Sum <- as.data.frame(table(defer_sc$ID))
table(Sum$Freq)


#--------------------------#
#### 4.3 Run the models ####
#--------------------------#

## run the binomial GLM for deferal without sub-population
def_model <- glmmTMB(attempt ~  Comp1 + Gr10precip + Green_centre + year + (1|ID),
                      data = defer_sc,
                      family = binomial(link = "logit"))
## examine model
summary(def_model)
drop1(def_model, test = "Chi") # Liklihood ratio test, no terms significant
confint(def_model)


## run the binomial GLM for defer with only sub-popualtion
def_comp <- glmmTMB(attempt ~  sub_pop + year + (1|ID),
                      data = defer_sc,
                      family = binomial(link = "logit"))

## examine model
summary(def_comp)
drop1(def_comp, test = "Chi") # Liklihood ratio test, no terms significant
confint(def_comp)


## test model assumptions using DHARMa
## simulate randomized quantile residuals
simulationOutput <- simulateResiduals(fittedModel = def_model_ran, plot = T) # not violated

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

## perform model averaging
# select models within 5 AICc points of top model
msdef_sub <- subset(msdef, !nested(.), recalc.weights=T)
msdef_sub <- subset(msdef_sub, delta <= 6, recalc.weights=T)
msdef_sub # only retained intercept only model



#--------------------------------#
#### 4.4 Plot model estimates ####
#--------------------------------#

## get the estimates and confidence intervals from the models
Ests <- as.data.frame(confint(def_model)) %>% 
        slice(., 1:(n()-1))
rownames(Ests) <- c("Intercept", "Arrival Climate PCA", "Arrival Precipitation", "Arrival Date",
                     "Year [2019]", "Year [2020]", "Year [2021]")

## make a forest plot fo the model estimates
Ests <- Ests %>% 
  map_df(rev) %>% 
  mutate(Param = rev(rownames(Ests)),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))

Ests$Param <- factor(Ests$Param, levels = Ests$Param[order(1:9)])

Def <- ggplot(Ests) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
  geom_errorbarh(aes(y=Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, size =0.5) +
  geom_point(aes(y=Param, x= Estimate, color = Sig), size = 2.5) +
  theme_bw() +
  ylab("") +
  xlim(-10.5, 11.4) +
  scale_color_manual(values=c("#B2BABB", "#E74C3C")) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title=element_text(size=16), 
        legend.title=element_text(size=14),
        axis.text=element_text(size=14), 
        legend.text=element_text(size=12),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

Def

## save the plot as a png
ggsave("Paper Plots/Supp Fig 8- Breeding Deferal forest plot.png",
       width = 22, height = 20, units = "cm")





#-------------------------------------------------------#
#### 5. Binomal GLM for succes/fail of attempted breeding
#-------------------------------------------------------#

#-------------------------------------#
#### 5.1 Prepare data for modeling ####
#-------------------------------------#

## filter out birds that did not attmept to breed
# breeders <- dplyr::filter(Inc_acc, is.na(attempt_start) == F)
breeders <- Inc_acc

## set response as correct class 
breeders$success24 <- ifelse(is.na(breeders$success24 ==T), 0, breeders$success24)
breeders$success24 <- as.factor(breeders$success24)
breeders$year <- as.factor(breeders$year)

## scale the explanatories
breed_sc <- breeders
breed_sc <- breeders %>% 
            mutate(Gr10precip = (Gr10precip-mean(Gr10precip, na.rm=T))/sd(Gr10precip, na.rm=T),
                   staging_length = (staging_length-mean(staging_length, na.rm=T))/sd(staging_length, na.rm=T),
                   Green_centre = (Green_centre-mean(Green_centre, na.rm=T))/sd(Green_centre, na.rm=T),
                   laying_centre = (laying_centre-mean(laying_centre, na.rm=T))/sd(laying_centre, na.rm=T),
                   breeding_lat = (breeding_lat-mean(breeding_lat, na.rm=T))/sd(breeding_lat, na.rm=T),
                   Comp1 = (Comp1-mean(Comp1, na.rm=T))/sd(Comp1, na.rm=T))

## check for correlations between possible explantory variables
test2 <- subset(breed_sc, select = c("Gr10precip", "breeding_lat", "staging_length", 
                                     "Green_centre", "laying_centre", "Comp1", "Comp2"))
correlation_matrix2 <- cor(test2) # Note: breedng lat has moderate negative correlation with two PCA axis

## assign sub-population
breed_sc$sub_pop <- ifelse(breed_sc$Ringing.location == "WEXF" | breed_sc$Ringing.location == "SESK", "HVAN", "LOWLANDS")

## Check if other variables differ between levels of sub-population
biserial.cor(breed_sc$Comp1, breed_sc$sub_pop)
biserial.cor(breed_sc$Gr10precip, breed_sc$sub_pop)
biserial.cor(breed_sc$Green_centre, breed_sc$sub_pop)



#-----------------------------------#
#### 5.2 Run model and do checks ####
#-----------------------------------#

## Model incubation success for both sub-populations together
breed_mod <- glmmTMB(success24 ~ Comp1 + Gr10precip + Green_centre + year + (1|ID),
                   data = breed_sc,
                   family = binomial(link = "logit"))
summary(breed_mod)
drop1(breed_mod, test = "Chi")
confint(breed_mod)

## Comparison model between sub-populations
breed_comp <- glmmTMB(success24 ~ sub_pop + year + (1|ID),
                       data = breed_sc,
                       family = binomial(link = "logit"))
summary(breed_comp)
drop1(breed_comp, test = "Chi")
confint(breed_comp)



## check model assumptions with DHARMa
## simulate randomized quantile residuals
simulationOutput24 <- simulateResiduals(fittedModel = breed_mod24_ran, plot = T) 

## plot the residuals against other predictors
par(mfrow=c(4,2))
plotResiduals(simulationOutput24, form = breed_sc$Gr10precip) # not violated
plotResiduals(simulationOutput24, form = breed_sc$breeding_lat) # not violated
plotResiduals(simulationOutput24, form = breed_sc$Green_centre) # not violated shows slight quadratic shape in plot but probaly fine to leave
plotResiduals(simulationOutput24, form = breed_sc$laying_centre) # not violated
plotResiduals(simulationOutput24, form = breed_sc$mig_duration) # not violated
plotResiduals(simulationOutput24, form = breed_sc$staging_length) 
plotResiduals(simulationOutput24, form = breed_sc$Comp2)

## test for outliers, dispersion and confomrity to expected distribution
testResiduals(simulationOutput24) # not violated




#--------------------------------#
#### 5.3 Plot model estimates ####
#--------------------------------#


## get the estimates and confidence intervals from the models
Ests2 <- as.data.frame(confint(breed_mod)) %>% 
         slice(., 1:(n()-1))
rownames(Ests2) <- c("Intercept", "Arrival Climate PCA", "Arrival Precipitation", "Arrival Date",
                     "Year [2019]", "Year [2020]", "Year [2021]")

## make a forest plot fo the model estimates
Ests2 <- Ests2 %>% 
  map_df(rev) %>% 
  mutate(Param = rev(rownames(Ests2)),
         Sig = ifelse( (`2.5 %` >0 & `97.5 %` >0) | (`2.5 %` <0 & `97.5 %` <0), "red", "grey"))

Ests2$Param <- factor(Ests2$Param, levels = Ests2$Param[order(1:9)])

Succ <- ggplot(Ests2) +
          geom_vline(xintercept = 0, linetype = "dashed", alpha =0.5) +
          geom_errorbarh(aes(y=Param, xmin= `2.5 %`, xmax=`97.5 %`), height = 0.2, size =0.5) +
          geom_point(aes(y=Param, x= Estimate, color = Sig), size = 2.5) +
          theme_bw() +
          ylab("") +
          xlim(-4.1, 4) +
          scale_color_manual(values=c("#B2BABB", "#E74C3C")) +
          theme(panel.grid.minor.y = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_blank(),
                axis.title=element_text(size=16), 
                legend.title=element_text(size=14),
                axis.text=element_text(size=14), 
                legend.text=element_text(size=12),
                panel.grid.minor.x = element_blank(),
                legend.position = "none")

Succ


## save the plot as a png
ggsave("Paper Plots/Supp Fig 9- Breeding success forest plot.png",
       width = 22, height = 20, units = "cm")





