## Luke Ozsanlav-Harris

## Summarize the temp and precipitation annotation of the Greenland nest sites
## Create average temp, summed precipitation for a pre breeding period, post breeding period and full period

## packages required
pacman::p_load(ggplot2, data.table, lubridate, zoo, MuMIn, glmmTMB, DHARMa,
               nlme, lmtest, effects, performance, ggpubr, patchwork, png)
library(dplyr)
library(purrr)




#-------------------------------------------#
## 1. Calculate breeding Phenology dates ####
#-------------------------------------------#

## read in incubation data
Inc <- fread("Outputs/Incubation_attempt_lengths_new.csv")

## calculate median arrival date
arrival <- median(yday(Inc$Green_arrive))

## calculate median laying date
laying <- median(yday(Inc$attempt_start[is.na(Inc$attempt_start) == F]))

## calculate median hatching date (use 24 days as cut off for success)
hatch <- laying + 25





#------------------------------------------#
#### 2. Read in and sort annotated data ####
#------------------------------------------#

## Read in the climate data
clim <- readRDS("Data/All_sampling_sites_annotated.RDS")

## shorten the column names
colnames(clim)
setnames(clim, old = c("ECMWF Interim Full Daily SFC Temperature (2 m above Ground)", "ECMWF Interim Full Daily SFC-FC Total Precipitation"), 
         new = c("temp", "precip"))

## set precip values below 0.00005 to zero
clim$precip <- ifelse(clim$precip < 0.00005, 0, clim$precip)

## Convert units of climate data
clim$temp <- clim$temp - 273.15 # convert to degrees centigrade
clim$precip <- clim$precip * 1000 # was measured in m so divide by 1000 to get to mm






#-------------------------------------------------------#
#### 3. Summarize temp and precip by breeding period ####
#-------------------------------------------------------#

## Add year and month columns
clim <- clim %>% mutate(year = year(clim$timestamp),
                        month = month (clim$timestamp),
                        yday = yday(clim$timestamp))


## Add column to discern two different periods of breedin
clim$t_period <- ifelse(clim$yday < hatch, "pre_hatch", "post_hatch")


## Aggregate the climate by day, year and sampling location
Clim_day <- clim %>% 
            group_by(Ringing_loc, t_period, year, yday, Source_ID) %>% 
            summarise(ave_temp = mean(temp), 
                      tot_precip = mean(precip))


## Add column to indicate if temp is below or above zero
Clim_day$sub_zero <- ifelse(Clim_day$ave_temp < 0, 1, 0)

## Create summary by month and region
## create a load if indicator columns were values are only non-zero if they are in the correct time period
Clim_day$pre_hatch_temp <- ifelse(Clim_day$t_period== "pre_hatch", Clim_day$ave_temp, 0)
Clim_day$post_hatch_temp <- ifelse(Clim_day$t_period== "post_hatch", Clim_day$ave_temp, 0)

Clim_day$pre_hatch_precip <- ifelse(Clim_day$t_period== "pre_hatch", Clim_day$tot_precip, 0)
Clim_day$post_hatch_precip <- ifelse(Clim_day$t_period== "post_hatch", Clim_day$tot_precip, 0)#

Clim_day$pre_hatch_freeze <- ifelse(Clim_day$t_period== "pre_hatch", Clim_day$sub_zero, 0)
Clim_day$post_hatch_freeze <- ifelse(Clim_day$t_period== "post_hatch", Clim_day$sub_zero, 0)


## now summarize these by ringing site and year
Clim_period <- Clim_day %>% 
                 group_by(Ringing_loc, year) %>% 
                 summarise(n_samples = length(unique(Source_ID)),
                           pre_hatch_temp = mean(pre_hatch_temp),
                           post_hatch_temp = mean(post_hatch_temp),
                           pre_hatch_precip = sum(pre_hatch_precip)/n_samples,
                           post_hatch_precip = sum(post_hatch_precip)/n_samples,
                           pre_hatch_freeze = sum(pre_hatch_freeze)/n_samples,
                           post_hatch_freeze = sum(post_hatch_freeze)/n_samples) %>% 
                 ungroup()

## create a number of additonal combos that are just cimposites of the columns that laready exists
Clim_period$All_temp <- (Clim_period$pre_hatch_temp + Clim_period$post_hatch_temp)/2
Clim_period$All_precip <- (Clim_period$pre_hatch_precip + Clim_period$post_hatch_precip)
Clim_period$All_freeze <- (Clim_period$pre_hatch_freeze + Clim_period$post_hatch_freeze)






#------------------------------------------------------#
#### 4. Read in the productivity and snow melt data ####
#------------------------------------------------------#

## read in producivity data 
Productiv <- read.csv("Data/GWfG Islay Wexford young counts long format.csv")

## change column names so that data sets will join together
colnames(Productiv)
setnames(Productiv, old = c("data_year", "Ringing.location"), new = c("year", "Ringing_loc"))



## Join productivity and climate data
Clim_period2 <- inner_join(Clim_period, Productiv, by = c("year", "Ringing_loc"))
stopifnot(nrow(Clim_period)==nrow(Clim_period2)) # check the join worked properly






#------------------------------------------------------------------#
#### 5. Plot the relationships between climate and productivity ####
#------------------------------------------------------------------#

## Plots of all the data
Clim_period2$pop_trend <- ifelse(Clim_period2$year > 1999, "pop decline", "pop increase")
# p1 <- ggplot(Clim_period2, aes(x= All_temp, y= Perc_young, colour = pop_trend)) + geom_point() + facet_wrap(~Ringing_loc + pop_trend )
# p2 <- ggplot(Clim_period2, aes(x= All_precip, y= Perc_young, colour = pop_trend)) + geom_point() + facet_wrap(~Ringing_loc + pop_trend )
# p3 <- ggplot(Clim_period2, aes(x= All_freeze, y= Perc_young, colour = pop_trend)) + geom_point() + facet_wrap(~Ringing_loc + pop_trend )





#-------------------------------------------------------------------#
#### 6. Model the relationships between climate and productivity ####
#-------------------------------------------------------------------#


#---------------------------------------------------------#
#### 6.1 Organize data sets and test autocorrelation ####
#---------------------------------------------------------#

## create a column that is number of young per 1000
Clim_period2$Young <- Clim_period2$Perc_young*10

## create column for pop_trend 
Clim_period2$Trend <- ifelse(Clim_period2$year <= 1999, "Inc", "Dec")

## SPlit up the data into Islay and Wexford
Clim_Islay <- filter(Clim_period2, Ringing_loc == "ISLA")
Clim_Wexf <- filter(Clim_period2, Ringing_loc == "WEXF")

## test for auto correlation
modeltest <- lm(Young ~ All_temp + All_precip + Trend, data = Clim_Islay)
modeltest2 <- lm(Young ~ All_temp + All_precip + Trend, data = Clim_Wexf)
lmtest::dwtest(modeltest); lmtest::dwtest(modeltest2)




#------------------------------------------------------------#
#### 6.2 Check for temporal trends in the number of young ####
#------------------------------------------------------------#

## Is there autocorrelation in the NULL model and maximal model
lmtest::dwtest(lm(Young ~ 1, data = Clim_period2))
lmtest::dwtest(lm(Young ~ year*Ringing_loc, data = Clim_period2))

## Run model without autocorrelation and with autocorrelation
## Compare the different models with AICc to see what fits better
csY <- corARMA(c(0.3, -0.3), p = 2, q = 0, form = ~year|Ringing_loc) # define a second order autocorrelation
TrendyNo <- gls(Young ~ year*Ringing_loc, data = Clim_period2)
Trendy <- gls(Young ~ year*Ringing_loc, correlation = corAR1(form = ~year|Ringing_loc), data = Clim_period2)
Trendy2 <- gls(Young ~ year*Ringing_loc, correlation = csY, data = Clim_period2)

## AICc comparison, the autocorrelation does not make a difference
## Will use the first order autocorrelation as there is autocorrelation in the wexford time series
AICc(TrendyNo, Trendy, Trendy2)

## plot and summarize the output
top_mod_effects <- effects::predictorEffects(Trendy); plot(top_mod_effects)
summary(Trendy)
confint(Trendy)

## get confidence intervals for both of the slopes
Clim_period2$Ringing_loc <- as.factor(Clim_period2$Ringing_loc)
Clim_period2$Ringing_loc <- relevel(Clim_period2$Ringing_loc, "WEXF")
Trendy2 <- gls(Young ~ year*Ringing_loc, correlation = corAR1(form = ~year|Ringing_loc), data = Clim_period2)
confint(Trendy2)
confint(Trendy)
summary(Trendy2)
summary(Trendy)


## Create the plot for publication....
## use the effects package to extract the fit for the first variable
divisions <- 200
Trendy_effects <- predictorEffects(Trendy, focal.levels = divisions)
plot(Trendy_effects[1])
effectsTrendy <- Trendy_effects[1]
fitTrendy <- as.data.frame(cbind(effectsTrendy[["year"]][["fit"]], effectsTrendy[["year"]][["lower"]], 
                               effectsTrendy[["year"]][["upper"]], effectsTrendy[["year"]][["x"]][["Ringing_loc"]],
                               effectsTrendy[["year"]][["x"]][["year"]]))
## change the names to something meaningful
setnames(fitTrendy, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "Ringing_loc", "year"))
fitTrendy$Ringing_loc <- ifelse(fitTrendy$Ringing_loc == 1, "Wexford (Ire)", "Islay (Scot)")


## Now plot using ggplot
setnames(fitTrendy,old = "fit", new = "Young")
Clim_period3 <- Clim_period2
Clim_period3$Ringing_loc <- ifelse(Clim_period3$Ringing_loc == "WEXF", "Wexford (Ire)", "Islay (Scot)")

ggplot(mapping=aes(x= year, y = Young, group = Ringing_loc, colour = Ringing_loc)) + 
  geom_point(data = Clim_period3, mapping = aes(x = year, y = Young, colour = Ringing_loc), size = 1.75, alpha = 0.5) +
  geom_ribbon(data = fitTrendy, mapping =aes(x=year, ymin = lower, ymax = upper, group = Ringing_loc), 
              alpha = 0.2, colour = NA, fill = "grey")+
  geom_line(data=fitTrendy, size = 1.25)  +
  xlab("Year") + ylab("Number of juveniles per 1000") +
  labs(colour="Sub-population") +
  theme_bw() +
  scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())

# ggsave("Paper Plots/Figure 6- Temporal trends in producitvity.png",
#        width = 25, height = 19, units = "cm")




#--------------------------------#
#### 6.3 Run Models for Islay ####
#--------------------------------#

## check for correlation between predictors
cor(Clim_Islay$pre_hatch_freeze, Clim_Islay$pre_hatch_precip)
cor(Clim_Islay$post_hatch_freeze, Clim_Islay$post_hatch_precip)
cor(Clim_Islay$All_freeze, Clim_Islay$All_precip)

## Check if we need to account for autocorrelation
lmtest::dwtest(lm(Young ~ 1, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ Trend*All_freeze + All_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ post_hatch_freeze*Trend + post_hatch_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ All_freeze*Trend + post_hatch_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ All_freeze*Trend + pre_hatch_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ post_hatch_freeze*Trend + All_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ pre_hatch_freeze*Trend + All_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ post_hatch_freeze*Trend + pre_hatch_precip*Trend, data = Clim_Islay))
lmtest::dwtest(lm(Young ~ pre_hatch_freeze*Trend + post_hatch_precip*Trend, data = Clim_Islay))


## Models for Islay onlys
Null <- gls(Young ~ 1, correlation = corAR1(), data = Clim_Islay, method = "ML")
IAll <- gls(Young ~ All_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")
IPre <- gls(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")
IPost <- gls(Young ~ post_hatch_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")

IAllPost <- gls(Young ~ All_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")
IAllPre <- gls(Young ~ All_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")

IPostAll <- gls(Young ~ post_hatch_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")
IPreAll <- gls(Young ~ pre_hatch_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")

IPostPre <- gls(Young ~ post_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")
IPrePost <- gls(Young ~ pre_hatch_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")


## compare the AICc of all the models
AICc(Null, IAll, IPre, IPost, IAllPost, IAllPre, IPostAll, IPreAll, IPostPre, IPrePost)
model_performance(IAll); model_performance(IPre); model_performance(IPost)


## For the best model, check the autocorrelation structure
csAll <- corARMA(c(0.3, -0.3), p = 2, q = 0, form = ~year) # define a second order autocorrelation
IAllNo <- gls(Young ~ All_freeze*Trend + All_precip*Trend,  data = Clim_Islay, method = "ML")
IAll1 <- gls(Young ~ All_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Islay, method = "ML")
IAll2 <- gls(Young ~ All_freeze*Trend + All_precip*Trend, csAll, data = Clim_Islay, method = "ML")
AICc(IAllNo, IAll1, IAll2)



## Make plots of the autocorrelations
par(mfrow=c(1,2))
E <- residuals(IAllNo, type = "normalized")
plot(acf(E, plot = F), main="Islay- No Autocorrelation Term")

E2 <- residuals(IAll1, type = "normalized")
plot(acf(E2, plot = F), main="Islay- First Order Autocorrelation Term")

## save autocorrelation plot
png(filename = "Paper Plots/Supp Fig 2- ACF plots Islay.png", res = 300,
    width = 650, height = 480)



## get the p-values for the top model and plot the effects
summary(IAll) # model summary
drop1(IAll, test = "Chisq") # liklihood ratio test
confint(IAll) # CIs
Itop_mod_effects <- effects::predictorEffects(IAll) # model effects
plot(Itop_mod_effects)
emmeans::emmeans(IAll, specs = ~Trend*All_freeze)

## get the other slope estimates
Clim_Islay$Trend <- relevel(as.factor(Clim_Islay$Trend), "Inc")
IAll2 <- gls(Young ~ Trend*All_freeze + All_precip*Trend, correlation = corAR1(), data = Clim_Islay, method = "ML")
summary(IAll2);confint(IAll2)

## Look for trends in the variable of interest
#Trend1 <- gls(All_freeze ~ year, correlation = corAR1(), data = Clim_Islay)
#top_mod_effects <- effects::predictorEffects(Trend1); plot(top_mod_effects)
#summary(Trend1)

## Use dredge on the top model
options(na.action = "na.fail") 
Ims2 <- MuMIn::dredge(IAll, trace = 2, rank = "AICc")
Ims2_sub <- subset(Ims2, !nested(.), recalc.weights=T) ## this bit of code stops nested models being included in the subst of models for av
Ims2_sub <- subset(Ims2_sub, delta < 6, recalc.weights=T)




#---------------------------------#
#### 6.4 Plot Models for Islay ####
#---------------------------------#

## Run the top model from MuMin and then plot it
ITop <- gls(Young ~ All_freeze*Trend + All_precip*Trend, correlation = corAR1(), data = Clim_Islay, method = "ML")
confint(ITop)


## use the effects package to extract the fit for the first variable
divisions <- 200
ITop_mod_effects <- predictorEffects(ITop, focal.levels = divisions)
plot(ITop_mod_effects[1])
effectsITop <- ITop_mod_effects[1]
fitITop <- as.data.frame(cbind(effectsITop[["All_freeze"]][["fit"]], effectsITop[["All_freeze"]][["lower"]], 
                               effectsITop[["All_freeze"]][["upper"]], effectsITop[["All_freeze"]][["x"]][["Trend"]],
                               effectsITop[["All_freeze"]][["x"]][["All_freeze"]]))
## change the names to something meaningful
setnames(fitITop, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "Trend", "All_freeze"))
fitITop$Trend <- ifelse(fitITop$Trend == 1, "Increasing", "Decreasing")


## Now plot using ggplot
setnames(fitITop,old = "fit", new = "Young")
Clim_Islay2 <- Clim_Islay
Clim_Islay2$Trend <- ifelse(Clim_Islay2$Trend == "Dec", "Decreasing", "Increasing")


I1 <- ggplot(mapping=aes(x= All_freeze, y = Young, group = Trend, colour = Trend)) + 
      geom_ribbon(data = fitITop, mapping =aes(x=All_freeze, ymin = lower, ymax = upper, group = Trend), 
                  alpha = 0.2, colour = NA, fill = "grey")+
      geom_line(data=fitITop, size = 1.25)  +
      geom_point(data = Clim_Islay2, alpha = 0.5) +
      xlab("Days below 0°C (April 1st–August 31st)") + ylab("Number of juveniles per 1000") +
      labs(colour="Global trend:") +
      theme_bw() +
      scale_colour_manual(values=c("#E1BE6A", "#40B0A6")) +
      theme(panel.grid.minor.y = element_blank(),
            axis.title=element_text(size=18),
            axis.text=element_text(size=14),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14),
            panel.grid.minor.x = element_blank(), 
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.spacing.x = unit(0.4, "cm")) + 
     annotate(geom="text", x=72, y=355, label="Islay", color="#0072B2", size =8)




## use the effects package to extract the fit for the first variable
divisions <- 200
ITop_mod_effects2 <- predictorEffects(ITop, focal.levels = divisions)
plot(ITop_mod_effects2[3])
effectsITop2 <- ITop_mod_effects2[3]
fitITop2 <- as.data.frame(cbind(effectsITop2[["All_precip"]][["fit"]], effectsITop2[["All_precip"]][["lower"]], 
                               effectsITop2[["All_precip"]][["upper"]], effectsITop2[["All_precip"]][["x"]][["Trend"]],
                               effectsITop2[["All_precip"]][["x"]][["All_precip"]]))
## change the names to something meaningful
setnames(fitITop2, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "Trend", "All_precip"))
fitITop2$Trend <- ifelse(fitITop2$Trend == 1, "Increasing", "Decreasing")


## Now plot using ggplot
setnames(fitITop2,old = "fit", new = "Young")
Clim_Islay2 <- Clim_Islay
Clim_Islay2$Trend <- ifelse(Clim_Islay2$Trend == "Dec", "Decreasing", "Increasing")

I2 <- ggplot(mapping=aes(x= All_precip, y = Young, group = Trend, colour = Trend)) + 
      geom_ribbon(data = fitITop2, mapping =aes(x=All_precip, ymin = lower, ymax = upper, group = Trend), 
                  alpha = 0.2, colour = NA, fill = "grey")+
      geom_line(data=fitITop2, size = 1.25)  +
      geom_point(data = Clim_Islay2, alpha = 0.5) +
      xlab("Total precipitation/mm (April 1st–August 31st)") + ylab("Number of juveniles per 1000") +
      labs(colour="Global trend:") +
      theme_bw() +
      scale_colour_manual(values=c("#E1BE6A", "#40B0A6")) +
      theme(panel.grid.minor.y = element_blank(),
            axis.title=element_text(size=18),
            axis.text=element_text(size=14),
            legend.text=element_text(size=14),
            legend.title=element_text(size=14),
            panel.grid.minor.x = element_blank(), 
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.spacing.x = unit(0.4, "cm")) +
      annotate(geom="text", x=84, y=285, label="Islay", color="#0072B2", size =8)




#----------------------------------#
#### 6.5 Run models for Wexford ####
#----------------------------------#

## check for correlations betwen variables
cor(Clim_Wexf$pre_hatch_freeze, Clim_Wexf$pre_hatch_precip)
cor(Clim_Wexf$post_hatch_freeze, Clim_Wexf$post_hatch_precip)
cor(Clim_Wexf$All_freeze, Clim_Wexf$All_precip)
cor(Clim_Wexf$All_freeze, Clim_Wexf$pre_hatch_precip)
cor(Clim_Wexf$All_freeze, Clim_Wexf$post_hatch_precip)
cor(Clim_Wexf$All_precip, Clim_Wexf$pre_hatch_freeze)

## Check that what the model without correlations is like
lmtest::dwtest(lm(Young ~ 1, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ Trend*All_freeze + All_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ post_hatch_freeze*Trend + post_hatch_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ All_freeze*Trend + post_hatch_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ All_freeze*Trend + pre_hatch_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ post_hatch_freeze*Trend + All_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ pre_hatch_freeze*Trend + All_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ post_hatch_freeze*Trend + pre_hatch_precip*Trend, data = Clim_Wexf))
lmtest::dwtest(lm(Young ~ pre_hatch_freeze*Trend + post_hatch_precip*Trend, data = Clim_Wexf))


## Models for Islay onlys
Null <- gls(Young ~ 1, correlation = corAR1(), data = Clim_Wexf, method = "ML")
WAll <- gls(Young ~ All_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")
WPre <- gls(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")
WPost <- gls(Young ~ post_hatch_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")

WAllPost <- gls(Young ~ All_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")
WAllPre <- gls(Young ~ All_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")

WPostAll <- gls(Young ~ post_hatch_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")
WPreAll <- gls(Young ~ pre_hatch_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")

WPostPre <- gls(Young ~ post_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")
WPrePost <- gls(Young ~ pre_hatch_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")


## compare the AWCc of all the models
AICc(Null, WAll, WPre, WPost, WAllPost, WAllPre, WPostAll, WPreAll, WPostPre, WPrePost)
model_performance(WPre)

## For the best model, check the autocorrelation structure
csAll <- corARMA(c(0.3, -0.3), p = 2, q = 0, form = ~year) # define a second order autocorrelation
WAllNo <- gls(Young ~ All_freeze*Trend + All_precip*Trend,  data = Clim_Wexf, method = "ML")
WAll1 <- gls(Young ~ All_freeze*Trend + All_precip*Trend, correlation = corAR1(form = ~year), data = Clim_Wexf, method = "ML")
WAll2 <- gls(Young ~ All_freeze*Trend + All_precip*Trend, csAll, data = Clim_Wexf, method = "ML")
AICc(WAllNo, WAll1, WAll2)

## Make plots of the autocorrelations
par(mfrow=c(1,2))
E <- residuals(WAllNo, type = "normalized")
plot(acf(E, plot = F), main="Wexford- No Autocorrelation Term")

E2 <- residuals(WAll1, type = "normalized")
plot(acf(E2, plot = F), main="Wexford- First Order Autocorrelation Term")

## save autocorrelation plot
png(filename = "Paper Plots/Supp Fig 3- ACF plots Wexf.png", res = 300,
    width = 650, height = 480)







## get the p-values for the top model and plot the effects
summary(WPre); confint(WPre)
drop1(WAll, test = "Chisq")
Wtop_mod_effects <- effects::predictorEffects(WAll)
plot(Wtop_mod_effects)

## get the other slope estimates
Clim_Wexf$Trend <- relevel(as.factor(Clim_Wexf$Trend), "Inc")
Pre2 <- gls(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(), data = Clim_Wexf, method = "ML")
summary(Pre2);confint(Pre2)


## Look for trends in the variable of interest
# Trend1 <- gls(All_freeze ~ year, correlation = corAR1(), data = Clim_Islay)
# top_mod_effects <- effects::predictorEffects(Trend1); plot(top_mod_effects)
# summary(Trend1)

## Use dredge on the top model
options(na.action = "na.fail") 
Wms2 <- MuMIn::dredge(WAll_ML, trace = 2, rank = "AICc")
Wms2_sub <- subset(Wms2, !nested(.), recalc.weights=T) ## this bit of code stops nested models being included in the subst of models for av
Wms2_sub <- subset(Wms2_sub, delta < 6, recalc.weights=T)




#-----------------------------------#
#### 6.6 Plot models for Wexford ####
#-----------------------------------#

## Run the top model from MuMin and then plot it
WTop <- gls(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(), data = Clim_Wexf, method = "ML")
confint(WTop)

## use the effects package to extract the fit for the first variable
divisions <- 200
Trendy_effects <- predictorEffects(WTop, focal.levels = divisions)
plot(Trendy_effects[1])
effectsWTop <- Trendy_effects[1]
fitWTop <- as.data.frame(cbind(effectsWTop[["pre_hatch_freeze"]][["fit"]], effectsWTop[["pre_hatch_freeze"]][["lower"]], 
                               effectsWTop[["pre_hatch_freeze"]][["upper"]], effectsWTop[["pre_hatch_freeze"]][["x"]][["Trend"]],
                               effectsWTop[["pre_hatch_freeze"]][["x"]][["pre_hatch_freeze"]]))
## change the names to something meaningful
setnames(fitWTop, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "Trend", "pre_hatch_freeze"))
fitWTop$Trend <- ifelse(fitWTop$Trend == 1, "Increasing", "Decreasing")


## Now plot using ggplot
setnames(fitWTop,old = "fit", new = "Young")
Clim_Wexf2 <- Clim_Wexf
Clim_Wexf2$Trend <- ifelse(Clim_Wexf2$Trend == "Dec", "Decreasing", "Increasing")

W1<-ggplot(mapping=aes(x= pre_hatch_freeze, y = Young, group = Trend, colour = Trend)) + 
    geom_ribbon(data = fitWTop, mapping =aes(x=pre_hatch_freeze, ymin = lower, ymax = upper, group = Trend), 
                alpha = 0.2, colour = NA, fill = "grey")+
    geom_line(data=fitWTop, size = 1.25)  +
    geom_point(data = Clim_Wexf2, alpha = 0.5) +
    xlab("Days below 0°C (April 1st–June 20th)") + ylab("Number of juveniles per 1000") +
    labs(colour="Global trend") +
    theme_bw() +
    scale_colour_manual(values=c("#E1BE6A", "#40B0A6")) +
    theme(panel.grid.minor.y = element_blank(),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14),
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.spacing.x = unit(0.4, "cm"))+ 
  annotate(geom="text", x=70, y=340, label="Wexford", color="#D55E00", size =8)





## use the effects package to extract the fit for the second variable
divisions <- 200
Trendy_effects2 <- predictorEffects(WTop, focal.levels = divisions)
plot(Trendy_effects2[3])
effectsWTop2 <- Trendy_effects2[3]
fitWTop2 <- as.data.frame(cbind(effectsWTop2[["pre_hatch_precip"]][["fit"]], effectsWTop2[["pre_hatch_precip"]][["lower"]], 
                               effectsWTop2[["pre_hatch_precip"]][["upper"]], effectsWTop2[["pre_hatch_precip"]][["x"]][["Trend"]],
                               effectsWTop2[["pre_hatch_precip"]][["x"]][["pre_hatch_precip"]]))
## change the names to something meaningful
setnames(fitWTop2, old = c("V1", "V2", "V3", "V4", "V5"), new = c("fit", "lower", "upper", "Trend", "pre_hatch_precip"))
fitWTop2$Trend <- ifelse(fitWTop2$Trend == 1, "Increasing", "Decreasing")


## Now plot using ggplot
setnames(fitWTop2,old = "fit", new = "Young")
Clim_Wexf2 <- Clim_Wexf
Clim_Wexf2$Trend <- ifelse(Clim_Wexf2$Trend == "Dec", "Decreasing", "Increasing")

W2<-ggplot(mapping=aes(x= pre_hatch_precip, y = Young, group = Trend, colour = Trend)) + 
    geom_ribbon(data = fitWTop2, mapping =aes(x=pre_hatch_precip, ymin = lower, ymax = upper, group = Trend), 
                alpha = 0.2, colour = NA, fill = "grey")+
    geom_line(data=fitWTop2, size = 1.25)  +
    geom_point(data = Clim_Wexf2, alpha = 0.5) +
    xlab("Total Precipitation/mm (April 1st–June 20th)") + ylab("Number of juveniles per 1000") +
    labs(colour="Global trend:") +
    theme_bw() +
    scale_colour_manual(values=c("#E1BE6A", "#40B0A6")) +
    theme(panel.grid.minor.y = element_blank(),
          axis.title=element_text(size=18),
          axis.text=element_text(size=14),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14),
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.spacing.x = unit(0.4, "cm")) +
  annotate(geom="text", x=48, y=340, label="Wexford", color="#D55E00", size =8)


## now arrange all of the plots onto a single plot
ggarrange(I1, W1, I2, W2, common.legend = T, nrow = 2, ncol = 2, legend = "top")

# ggsave("Paper Plots/Figure 7- Freezedays and Precip vs producitvity.png",
#        width = 32, height = 26, units = "cm")





#----------------------------------#
## 7. Temporal trends in climate ###
#----------------------------------#

#--------------------------------------#
#### 7.1 Plot trends in Temperature ####
#--------------------------------------#

## Trend in Islay Temperature
Trend1 <- gls(All_freeze ~ year, correlation = corAR1(), data = Clim_Islay)
summary(Trend1)
anova(Trend1)
confint(Trend1)

## extract the model fit
divisions <- 200
IslayTemp <- predictorEffects(Trend1, focal.levels = divisions)
plot(IslayTemp[1])
effectsIslay1 <- IslayTemp[1]
fitIslayT <- as.data.frame(cbind(effectsIslay1[["year"]][["fit"]], effectsIslay1[["year"]][["lower"]], 
                                 effectsIslay1[["year"]][["upper"]], effectsIslay1[["year"]][["x"]][["year"]]))

## change the names to something meaningful
setnames(fitIslayT, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "year"))

## make the plot
IT <- ggplot() + 
        geom_point(data = Clim_Islay, aes(x = year, y = All_freeze), shape = 1, colour = "grey", stroke = 1.5) +
        geom_ribbon(data = fitIslayT, mapping =aes(x=year, ymin = lower, ymax = upper), alpha = 0.2, colour = NA, fill = "grey")+
        geom_line(data=fitIslayT, aes(x= year, y = fit), size = 1.25, colour = "#0072B2")  +
        xlab("Year") + ylab("Days Sub 0°C (April 1st–August 31st)") +
        theme_bw() + ylim(28, 80) +
        #scale_colour_manual(values=c("#0072B2", "#D55E00")) +
        theme(panel.grid.minor.y = element_blank(),
              axis.title=element_text(size=18),
              axis.text=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14),
              panel.grid.minor.x = element_blank(), 
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
        annotate(geom="text", x= 2017.5, y= 79, label="Islay*", color="#0072B2", size =8)




## Trend in Wexford Temp
Trend2 <- gls(pre_hatch_freeze ~ year, correlation = corAR1(), data = Clim_Wexf)
summary(Trend2)
anova(Trend2)
confint(Trend2)

## extract the model fit
divisions <- 200
WexfTemp <- predictorEffects(Trend2, focal.levels = divisions)
plot(WexfTemp[1])
effectsWexf1 <- WexfTemp[1]
fitWexfT <- as.data.frame(cbind(effectsWexf1[["year"]][["fit"]], effectsWexf1[["year"]][["lower"]], 
                                 effectsWexf1[["year"]][["upper"]], effectsWexf1[["year"]][["x"]][["year"]]))

## change the names to something meaningful
setnames(fitWexfT, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "year"))

## make the plot
WT <- ggplot() + 
        geom_point(data = Clim_Wexf, aes(x = year, y = All_freeze), shape = 1, colour = "grey", stroke = 1.5) +
        geom_ribbon(data = fitWexfT, mapping =aes(x=year, ymin = lower, ymax = upper), alpha = 0.2, colour = NA, fill = "grey")+
        geom_line(data=fitWexfT, aes(x= year, y = fit), size = 1.25, colour = "#D55E00")  +
        xlab("Year") + ylab("Days Sub 0°C (April 1st–June 20th)") +
        theme_bw() + ylim(28, 80) +
        #scale_colour_manual(values=c("#0072B2", "#D55E00")) +
        theme(panel.grid.minor.y = element_blank(),
              axis.title=element_text(size=18),
              axis.text=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14),
              panel.grid.minor.x = element_blank(), 
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      annotate(geom="text", x= 2015.5, y= 79, label="Wexford*", color="#D55E00", size =8)




#----------------------------------------#
#### 7.2 Plot trends in Precipitation ####
#----------------------------------------#

## Trend in Islay Precipitation
Trend1_2 <- gls(All_precip ~ year, correlation = corAR1(), data = Clim_Islay)
summary(Trend1_2)
anova(Trend1_2)
confint(Trend1_2)

## extract the model fit
divisions <- 200
IslayPrecip <- predictorEffects(Trend1_2, focal.levels = divisions)
plot(IslayPrecip[1])
effectsIslay2 <- IslayPrecip[1]
fitIslayP <- as.data.frame(cbind(effectsIslay2[["year"]][["fit"]], effectsIslay2[["year"]][["lower"]], 
                                 effectsIslay2[["year"]][["upper"]], effectsIslay2[["year"]][["x"]][["year"]]))

## change the names to something meaningful
setnames(fitIslayP, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "year"))

## make the plot
IP <- ggplot() + 
  geom_point(data = Clim_Islay, aes(x = year, y = All_precip), shape = 1, colour = "grey", stroke = 1.5) +
  geom_ribbon(data = fitIslayP, mapping =aes(x=year, ymin = lower, ymax = upper), alpha = 0.2, colour = NA, fill = "grey")+
  geom_line(data=fitIslayP, aes(x= year, y = fit), size = 1.25, colour = "#0072B2")  +
  xlab("Year") + ylab("Total Precip/mm (April 1st–August 31st)") +
  theme_bw() + ylim(10, 90) +
  #scale_colour_manual(values=c("#0072B2", "#D55E00")) +
  theme(panel.grid.minor.y = element_blank(),
        axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  annotate(geom="text", x= 2017.5, y= 89, label="Islay*", color="#0072B2", size =8)




## Trend in Wexford Precipitation
Trend2_2 <- gls(pre_hatch_precip ~ year, correlation = corAR1(), data = Clim_Wexf)
summary(Trend2_2)
anova(Trend2_2)
confint(Trend2_2)

## extract the model fit
divisions <- 200
WexfPrecip <- predictorEffects(Trend2_2, focal.levels = divisions)
plot(WexfPrecip[1])
effectsWexf2 <- WexfPrecip[1]
fitWexfP <- as.data.frame(cbind(effectsWexf2[["year"]][["fit"]], effectsWexf2[["year"]][["lower"]], 
                                effectsWexf2[["year"]][["upper"]], effectsWexf2[["year"]][["x"]][["year"]]))

## change the names to something meaningful
setnames(fitWexfP, old = c("V1", "V2", "V3", "V4"), new = c("fit", "lower", "upper", "year"))

## make the plot
WP <- ggplot() + 
        geom_point(data = Clim_Wexf, aes(x = year, y = pre_hatch_precip), shape = 1, colour = "grey", stroke = 1.5) +
        geom_ribbon(data = fitWexfP, mapping =aes(x=year, ymin = lower, ymax = upper), alpha = 0.2, colour = NA, fill = "grey")+
        geom_line(data=fitWexfP, aes(x= year, y = fit), size = 1.25, colour = "#D55E00")  +
        xlab("Year") + ylab("Total Precip/mm (April 1st–June 20th)") +
        theme_bw() + ylim(10, 90) +
        #scale_colour_manual(values=c("#0072B2", "#D55E00")) +
        theme(panel.grid.minor.y = element_blank(),
              axis.title=element_text(size=18),
              axis.text=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_text(size=14),
              panel.grid.minor.x = element_blank(), 
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
        annotate(geom="text", x= 2015.5, y= 89, label="Wexford", color="#D55E00", size =8)








#-----------------------------#
#### 7.3 Combine all plots ####
#-----------------------------#

## combine all of the plots
ggarrange(IT, WT, IP, WP, nrow=2, ncol = 2)


## save the plot
ggsave("Paper Plots/Supp Fig 5- Climatic trends.png",
       width = 36, height = 30, units = "cm")




