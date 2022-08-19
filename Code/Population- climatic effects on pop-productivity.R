## Luke Ozsanlav-Harris

## Summarize the temp and precipitation annotation of the Greenland nest sites
## Create average temp, summed precipitation for a pre breeding period, post breeding period and full period

## packages required
pacman::p_load(ggplot2, data.table, lubridate, zoo, MuMIn, glmmTMB, 
               nlme, lmtest, effects, performance, ggpubr)
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
clim <- fread("Data/All_sampling_sites_annotated.csv")

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
#### 6.1 Organize data sets and do data quality checks ####
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

## This is the model, accounts for autocorrelation
Trendy <- gls(Young ~ year*Ringing_loc, correlation = corAR1(form = ~year|Ringing_loc), data = Clim_period2)

## plot and summarize the output
top_mod_effects <- effects::predictorEffects(Trendy); plot(top_mod_effects)
summary(Trendy)
confint(Trendy)


## get confidecne intervals for both of the slopes
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
fitTrendy$Ringing_loc <- ifelse(fitTrendy$Ringing_loc == 1, "Islay (Scot)", "Wexford (Ire)")


## Now plot using ggplot
setnames(fitTrendy,old = "fit", new = "Young")
Clim_period3 <- Clim_period2
Clim_period3$Ringing_loc <- ifelse(Clim_period3$Ringing_loc == "WEXF", "Wexford (Ire)", "Islay (Scot)")

ggplot(mapping=aes(x= year, y = Young, group = Ringing_loc, colour = Ringing_loc)) + 
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

ggsave("Paper Plots/Figure 5- Temporal trends in producitvity.png",
       width = 25, height = 15, units = "cm")




#--------------------------------#
#### 6.3 Run Models for Islay ####
#--------------------------------#

## check for correlation between predictors
cor(Clim_Islay$pre_hatch_freeze, Clim_Islay$pre_hatch_precip)
cor(Clim_Islay$post_hatch_freeze, Clim_Islay$post_hatch_precip)
cor(Clim_Islay$All_freeze, Clim_Islay$All_precip)

## Models for Islay onlys
Null <- gls(Young ~ 1, correlation = corAR1(), data = Clim_Islay, method = "ML")
IAll <- gls(Young ~ Trend*All_freeze + All_precip*Trend, correlation = corAR1(), data = Clim_Islay, method = "ML")
IPre <- gls(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(), data = Clim_Islay, method = "ML")
IPost <- gls(Young ~ post_hatch_temp*Trend + post_hatch_precip*Trend, correlation = corAR1(), data = Clim_Islay, method = "ML")


## compare the AICc of all the models
AICc(Null, IAll, IPre, IPost)
model_performance(IAll); model_performance(IPre); model_performance(IPost)


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
      xlab("Days below 0°C (April 1st – August 31st)") + ylab("Number of juveniles per 1000; Islay") +
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
            legend.spacing.x = unit(0.4, "cm"))+ 
  annotate(geom="text", x=75, y=340, label="a)",color="black", size =8)


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
      xlab("Total precipitation/mm (April 1st – August 31st)") + ylab("Number of juveniles per 1000; Islay") +
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
            annotate(geom="text", x=90, y=300, label="b)",color="black", size =8)

## now arrange both of the plots onto a single plot
ggarrange(I1, I2, common.legend = T, nrow = 2)

ggsave("Paper Plots//Figure 6- Islay Freezedays and precip days vs producitvity.png", 
       width = 25, height = 30, units = "cm")



#----------------------------------#
#### 6.4 Run models for Wexford ####
#----------------------------------#

cor(Clim_Wexf$pre_hatch_freeze, Clim_Wexf$pre_hatch_precip)
cor(Clim_Wexf$post_hatch_freeze, Clim_Wexf$post_hatch_precip)
cor(Clim_Wexf$All_freeze, Clim_Wexf$All_precip)
cor(Clim_Wexf$All_freeze, Clim_Wexf$pre_hatch_precip)
cor(Clim_Wexf$All_freeze, Clim_Wexf$post_hatch_precip)
cor(Clim_Wexf$All_precip, Clim_Wexf$pre_hatch_freeze)

## Now do the same for the Wexford population
Null <- gls(Young ~ 1, correlation = corAR1(), data = Clim_Wexf)
WAll <- gls(Young ~ All_freeze*Trend + All_precip*Trend, correlation = corAR1(), data = Clim_Wexf, method = "ML")
WPre <- gls(Young ~ pre_hatch_freeze*Trend + pre_hatch_precip*Trend, correlation = corAR1(), data = Clim_Wexf, method = "ML")
WPost <- gls(Young ~ post_hatch_freeze*Trend + post_hatch_precip*Trend, correlation = corAR1(), data = Clim_Wexf, method = "ML")


## compare the AICc of all the models
AICc(Null, WAll, WPre, WPost, WAllPost)
summary(WPre)

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
#Trend1 <- gls(All_freeze ~ year, correlation = corAR1(), data = Clim_Islay)
#top_mod_effects <- effects::predictorEffects(Trend1); plot(top_mod_effects)
#summary(Trend1)

## Use dredge on the top model
options(na.action = "na.fail") 
Wms2 <- MuMIn::dredge(WAll_ML, trace = 2, rank = "AICc")
Wms2_sub <- subset(Wms2, !nested(.), recalc.weights=T) ## this bit of code stops nested models being included in the subst of models for av
Wms2_sub <- subset(Wms2_sub, delta < 6, recalc.weights=T)



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
    xlab("Days below 0°C (April 1st – June 20th)") + ylab("Number of juveniles per 1000; Wexford") +
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
  annotate(geom="text", x=73, y=340, label="a)",color="black", size =8)





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
    xlab("Total Precipitation/mm (April 1st – June 20th)") + ylab("Number of juveniles per 1000; Wexford") +
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
  annotate(geom="text", x=52, y=340, label="b)",color="black", size =8)


## now arrange both of the plots onto a single plot
ggarrange(W1, W2, common.legend = T, nrow = 2, legend = "top")

ggsave("Paper Plots//Figure 7- Wexford Freezedays and Precip vs producitvity.png", 
       width = 26, height = 30, units = "cm")



